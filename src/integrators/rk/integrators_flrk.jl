
using ForwardDiff

"Parameters for right-hand side function of formal Lagrangian Runge-Kutta methods."
mutable struct ParametersFLRK{DT,TT,VT,D,S} <: Parameters{DT,TT}
    v::VT
    Δt::TT

    a::Matrix{TT}
    â::Matrix{TT}
    c::Vector{TT}

    t::TT

    q::Vector{DT}
end

function ParametersFLRK(DT, D, v::VT, Δt::TT, tab) where {TT,VT}
    ParametersFLRK{DT,TT,VT,D,tab.s}(v, Δt, tab.a, tab.â, tab.c, 0, zeros(DT,D))
end

struct NonlinearFunctionCacheFLRK{DT}
    Q::Matrix{DT}
    V::Matrix{DT}
    Y::Matrix{DT}

    function NonlinearFunctionCacheFLRK{DT}(d, s) where {DT}

        # create internal stage vectors
        Q = zeros(DT,d,s)
        V = zeros(DT,d,s)
        Y = zeros(DT,d,s)

        new(Q, V, Y)
    end
end

@generated function compute_stages!(x::Vector{ST}, Q::Matrix{ST}, V::Matrix{ST}, Y::Matrix{ST},
                                    params::ParametersFLRK{DT,TT,VT,D,S}) where {ST,DT,TT,VT,D,S}

    tQ::Vector{ST} = zeros(ST,D)
    tV::Vector{ST} = zeros(ST,D)

    quote
        local tᵢ::TT

        @assert D == size(Q,1) == size(V,1) == size(Y,1)
        @assert S == size(Q,2) == size(V,2) == size(Y,2)

        # copy x to Y and compute Q = q + Δt Y
        for i in 1:size(Y,2)
            for k in 1:size(Y,1)
                Y[k,i] = x[D*(i-1)+k]
                Q[k,i] = params.q[k] + params.Δt * Y[k,i]
            end
        end

        # compute V = v(Q)
        for i in 1:S
            tᵢ = params.t + params.Δt * params.c[i]
            simd_copy_xy_first!($tQ, Q, i)
            params.v(tᵢ, $tQ, $tV)
            simd_copy_yx_first!($tV, V, i)
        end
    end
end

"Compute stages of formal Lagrangian Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersFLRK{DT,TT,VT,D,S}) where {ST,DT,TT,VT,D,S}

    cache = NonlinearFunctionCacheFLRK{ST}(D, S)

    quote
        compute_stages!(x, $cache.Q, $cache.V, $cache.Y, params)

        local y1::ST
        local y2::ST

        # compute b = - (Y-AV)
        for i in 1:S
            for k in 1:D
                y1 = 0
                y2 = 0
                for j in 1:S
                    y1 += params.a[i,j] * $cache.V[k,j]
                    y2 += params.â[i,j] * $cache.V[k,j]
                end
                b[D*(i-1)+k] = - $cache.Y[k,i] + (y1 + y2)
            end
        end
    end
end


"Formal Lagrangian Runge-Kutta integrator."
struct IntegratorFLRK{DT, TT, AT, FT, GT, VT, ΩT, dHT, SPT, ST, IT, N} <: Integrator{DT,TT}
    equation::VODE{DT,TT,AT,FT,GT,VT,ΩT,dHT,N}
    tableau::TableauFIRK{TT}
    Δt::TT

    params::SPT
    solver::ST
    iguess::InitialGuessODE{DT, TT, VT, IT}

    q::Vector{Vector{Double{DT}}}
    p::Vector{Vector{Double{DT}}}

    v::Vector{DT}
    y::Vector{DT}

    Q::Matrix{DT}
    V::Matrix{DT}
    P::Matrix{DT}
    F::Matrix{DT}
    G::Matrix{DT}
    ϑ::Matrix{DT}
    Y::Matrix{DT}
    Z::Matrix{DT}
    J::Array{DT,3}
end

function IntegratorFLRK(equation::VODE{DT,TT,AT,FT,GT,VT,ΩT,dHT,N}, tableau::TableauFIRK{TT}, Δt::TT) where {DT,TT,AT,FT,GT,VT,ΩT,dHT,N}
    D = equation.d
    M = equation.n
    S = tableau.q.s

    # create solution vector for internal stages / nonlinear solver
    x = zeros(DT, D*S)

    # create solution vectors
    q = Array{Vector{Double{DT}}}(M)
    p = Array{Vector{Double{DT}}}(M)
    for i in 1:M
        q[i] = zeros(Double{DT},D)
        p[i] = zeros(Double{DT},D)
    end

    # create velocity and update vector
    v = zeros(DT,D)
    y = zeros(DT,D)

    # create internal stage vectors
    Q = zeros(DT,D,S)
    V = zeros(DT,D,S)
    P = zeros(DT,D,S)
    F = zeros(DT,D,S)
    G = zeros(DT,D,S)
    ϑ = zeros(DT,D,S)
    Y = zeros(DT,D,S)
    Z = zeros(DT,D,S)
    J = zeros(DT,D,D,S)

    # create params
    params = ParametersFLRK(DT, D, equation.v, Δt, tableau.q)

    # create rhs function for nonlinear solver
    function_stages = (x,b) -> function_stages!(x, b, params)

    # create solver
    solver = get_config(:nls_solver)(x, function_stages)

    # create initial guess
    iguess = InitialGuessODE(get_config(:ig_interpolation), equation, Δt)

    # create integrator
    IntegratorFLRK{DT, TT, AT, FT, GT, VT, ΩT, dHT, typeof(params), typeof(solver), typeof(iguess.int), N}(
                                        equation, tableau, Δt, params, solver, iguess,
                                        q, p, v, y, Q, V, P, F, G, ϑ, Y, Z, J)
end


function initialize!(int::IntegratorFLRK{DT,TT}, sol::SolutionPODE, m::Int) where {DT,TT}
    @assert m ≥ 1
    @assert m ≤ sol.ni

    # copy initial conditions from solution
    get_initial_conditions!(sol, int.q[m], int.p[m], m)

    # initialise initial guess
    initialize!(int.iguess, m, sol.t[0], int.q[m])
end

function initial_guess!(int::IntegratorFLRK, m::Int)
    # compute initial guess for internal stages
    for i in 1:int.tableau.q.s
        evaluate!(int.iguess, m, int.y, int.v, int.tableau.q.c[i])
        for k in 1:int.equation.d
            int.V[k,i] = int.v[k]
        end
    end
    for i in 1:int.tableau.q.s
        for k in 1:int.equation.d
            int.solver.x[int.equation.d*(i-1)+k] = 0
            for j in 1:int.tableau.q.s
                int.solver.x[int.equation.d*(i-1)+k] += int.tableau.q.a[i,j] * int.V[k,j]
            end
        end
    end
end


"Integrate ODE with fully implicit Runge-Kutta integrator."
function integrate_step!(int::IntegratorFLRK{DT,TT}, sol::SolutionPODE{DT,TT}, m::Int, n::Int) where {DT,TT}
    @assert m ≥ 1
    @assert m ≤ sol.ni

    @assert n ≥ 1
    @assert n ≤ sol.ntime

    # set time for nonlinear solver
    int.params.t  = sol.t[0] + (n-1)*int.Δt
    int.params.q .= int.q[m]

    # compute initial guess
    initial_guess!(int, m)

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params, n)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params, n)

    # compute vector field at internal stages
    compute_stages!(int.solver.x, int.Q, int.V, int.Y, int.params)

    # compute final update for q
    update_solution!(int.q[m], int.V, int.tableau.q.b, int.tableau.q.b̂, int.Δt)

    # create temporary arrays
    tQ = zeros(DT, int.equation.d)
    tV = zeros(DT, int.equation.d)
    tϑ = zeros(DT, int.equation.d)
    tF = zeros(DT, int.equation.d)
    tJ = zeros(DT, int.equation.d, int.equation.d)

    te  = zeros(DT, int.equation.d)
    dQ1 = zeros(DT, int.equation.d)
    dQ2 = zeros(DT, int.equation.d)
    dV1 = zeros(DT, int.equation.d)
    dV2 = zeros(DT, int.equation.d)

    # println()

    # compute ϑ = α(Q)
    for i in 1:int.tableau.q.s
        tᵢ = int.params.t + int.Δt * int.tableau.q.c[i]
        simd_copy_xy_first!(tQ, int.Q, i)
        int.equation.α(tᵢ, tQ, tϑ)
        simd_copy_yx_first!(tϑ, int.ϑ, i)
    end

    # compute f_0(Q, V(Q)) = int.equation.f(t, Q, V, F)
    for i in 1:int.tableau.q.s
        tᵢ = int.params.t + int.Δt * int.tableau.q.c[i]
        simd_copy_xy_first!(tQ, int.Q, i)
        simd_copy_xy_first!(tV, int.V, i)
        int.equation.f(tᵢ, tQ, tV, tF)
        simd_copy_yx_first!(tF, int.F, i)
    end

    # compute Jacobian of v via ForwardDiff
    ϵ = sqrt(eps())
    for i in 1:int.tableau.q.s
        tᵢ = int.params.t + int.Δt * int.tableau.q.c[i]
        v_rev! = (y,x) -> int.equation.v(tᵢ,x,y)
        simd_copy_xy_first!(tQ, int.Q, i)
        ForwardDiff.jacobian!(tJ, v_rev!, tV, tQ)
        simd_copy_yx_first!(tJ, int.J, i)

        # println(tJ)
        # println(eig(tJ))
    end

    # set initial guess P = ϑ = α(Q)
    int.P .= int.ϑ

    # fixed-point iteration for P
    P  = zeros(int.P)
    δP = zeros(int.P)
    for r in 1:get_config(:nls_nmax)
        # copy P
        P .= int.P

        # contract J with (ϑ-p) and add to Z
        int.G .= 0
        for l in 1:int.tableau.q.s
            for i in 1:int.equation.d
                for j in 1:int.equation.d
                    int.G[i,l] += (int.ϑ[j,l] - int.P[j,l]) * int.J[j,i,l]
                end
            end
        end

        # set Z = A(F+G)
        int.Z .= 0
        for l in 1:int.tableau.q.s
            for i in 1:int.equation.d
                for k in 1:int.tableau.q.s
                    int.Z[i,k] += int.tableau.q.a[k,l] * (int.F[i,l] + int.G[i,l])
                end
            end
        end

        # update P
        for l in 1:int.tableau.q.s
            for i in 1:int.equation.d
                int.P[i,l] = int.p[m][i] + int.Δt * int.Z[i,l]
            end
        end

        δP .= (int.P .- P)

        # println(r, ", ", norm(δP) / norm(int.P))
        if norm(δP) / norm(int.P) ≤ get_config(:nls_stol)
            break
        end
    end

    # update G
    int.G .= 0
    for l in 1:int.tableau.q.s
        for i in 1:int.equation.d
            for j in 1:int.equation.d
                int.G[i,l] += (int.ϑ[j,l] - int.P[j,l]) * int.J[j,i,l]
            end
        end
    end

    # println(int.F)
    # println(int.G)

    # compute final update for p
    update_solution!(int.p[m], int.F, int.tableau.q.b, int.tableau.q.b̂, int.Δt)
    update_solution!(int.p[m], int.G, int.tableau.q.b, int.tableau.q.b̂, int.Δt)

    # println(int.p[m])

    # copy solution to initial guess
    update!(int.iguess, m, sol.t[0] + n*int.Δt, int.q[m])

    # take care of periodic solutions
    cut_periodic_solution!(int.q[m], int.equation.periodicity)

    # copy to solution
    copy_solution!(sol, int.q[m], int.p[m], n, m)
end
