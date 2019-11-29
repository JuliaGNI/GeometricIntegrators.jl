
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
    p::Vector{DT}
end

function ParametersFLRK(DT, D, v::VT, Δt::TT, tab) where {TT,VT}
    ParametersFLRK{DT,TT,VT,D,tab.s}(v, Δt, tab.a, tab.â, tab.c, 0, zeros(DT,D), zeros(DT,D))
end

struct NonlinearFunctionCacheFLRK{DT}
    Q::Vector{Vector{DT}}
    V::Vector{Vector{DT}}
    Y::Vector{Vector{DT}}

    function NonlinearFunctionCacheFLRK{DT}(d, s) where {DT}

        # create internal stage vectors
        Q = create_internal_stage_vector(DT, d, s)
        V = create_internal_stage_vector(DT, d, s)
        Y = create_internal_stage_vector(DT, d, s)

        new(Q, V, Y)
    end
end

function compute_stages!(x::Vector{ST}, Q::Vector{Vector{ST}}, V::Vector{Vector{ST}}, Y::Vector{Vector{ST}},
                                    params::ParametersFLRK{DT,TT,VT,D,S}) where {ST,DT,TT,VT,D,S}
    local tᵢ::TT

    @assert S == length(Q) == length(V) == length(Y)

    for i in 1:S
        @assert D == length(Q[i]) == length(V[i]) == length(Y[i])
        tᵢ = params.t + params.Δt * params.c[i]

        # copy x to Y
        for k in 1:D
            Y[i][k] = x[D*(i-1)+k]
        end

        # compute Q = q + Δt Y
        Q[i] .= params.q .+ params.Δt .* Y[i]

        # compute V = v(Q)
        params.v(tᵢ, Q[i], V[i])
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
                    y1 += params.a[i,j] * $cache.V[j][k]
                    y2 += params.â[i,j] * $cache.V[j][k]
                end
                b[D*(i-1)+k] = - $cache.Y[i][k] + (y1 + y2)
            end
        end
    end
end


"Formal Lagrangian Runge-Kutta integrator."
struct IntegratorFLRK{DT, TT, AT, FT, GT, VT, ΩT, dHT, SPT, ST, IT <: InitialGuessODE{DT,TT,VT}, N} <: DeterministicIntegrator{DT,TT}
    equation::VODE{DT,TT,AT,FT,GT,VT,ΩT,dHT,N}
    tableau::TableauFIRK{TT}
    Δt::TT

    params::SPT
    solver::ST
    iguess::IT

    q::Vector{Vector{TwicePrecision{DT}}}
    p::Vector{Vector{TwicePrecision{DT}}}

    v::Vector{DT}
    y::Vector{DT}

    Q::Vector{Vector{DT}}
    V::Vector{Vector{DT}}
    P::Vector{Vector{DT}}
    F::Vector{Vector{DT}}
    G::Vector{Vector{DT}}
    ϑ::Vector{Vector{DT}}
    Y::Vector{Vector{DT}}
    Z::Vector{Vector{DT}}

    J::Vector{Matrix{DT}}
    A::Array{DT,2}
end

function IntegratorFLRK(equation::VODE{DT,TT,AT,FT,GT,VT,ΩT,dHT,N}, tableau::TableauFIRK{TT}, Δt::TT) where {DT,TT,AT,FT,GT,VT,ΩT,dHT,N}
    D = equation.d
    M = equation.n
    S = tableau.q.s

    # create solution vector for internal stages / nonlinear solver
    x = zeros(DT, D*S)

    # create solution vectors
    q = create_solution_vector(DT, D, M)
    p = create_solution_vector(DT, D, M)

    # create velocity and update vector
    v = zeros(DT,D)
    y = zeros(DT,D)

    # create internal stage vectors
    Q = create_internal_stage_vector(DT, D, S)
    V = create_internal_stage_vector(DT, D, S)
    P = create_internal_stage_vector(DT, D, S)
    F = create_internal_stage_vector(DT, D, S)
    G = create_internal_stage_vector(DT, D, S)
    ϑ = create_internal_stage_vector(DT, D, S)
    Y = create_internal_stage_vector(DT, D, S)
    Z = create_internal_stage_vector(DT, D, S)

    J = Array{Matrix{DT}}(undef, S)
    # J = Vector{Matrix{DT}}(S)
    for i in 1:S
        J[i] = zeros(DT,D,D)
    end

    A = zeros(DT,D*S,D*S)


    # create params
    params = ParametersFLRK(DT, D, equation.v, Δt, tableau.q)

    # create rhs function for nonlinear solver
    function_stages = (x,b) -> function_stages!(x, b, params)

    # create solver
    solver = get_config(:nls_solver)(x, function_stages)

    # create initial guess
    iguess = InitialGuessODE(get_config(:ig_interpolation), equation, Δt)

    # create integrator
    IntegratorFLRK{DT, TT, AT, FT, GT, VT, ΩT, dHT, typeof(params), typeof(solver), typeof(iguess), N}(
                                        equation, tableau, Δt, params, solver, iguess,
                                        q, p, v, y, Q, V, P, F, G, ϑ, Y, Z, J, A)
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
            int.V[i][k] = int.v[k]
        end
    end
    for i in 1:int.tableau.q.s
        for k in 1:int.equation.d
            int.solver.x[int.equation.d*(i-1)+k] = 0
            for j in 1:int.tableau.q.s
                int.solver.x[int.equation.d*(i-1)+k] += int.tableau.q.a[i,j] * int.V[j][k]
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
    int.params.p .= int.p[m]

    # compute initial guess
    initial_guess!(int, m)

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params)

    # compute vector field at internal stages
    compute_stages!(int.solver.x, int.Q, int.V, int.Y, int.params)

    # compute final update for q
    update_solution!(int.q[m], int.V, int.tableau.q.b, int.tableau.q.b̂, int.Δt)

    # create temporary arrays
    tV = zeros(DT, int.equation.d)
    δP = zeros(DT, int.equation.d*int.tableau.q.s)

    # compute ϑ = α(Q), V(Q) = int.equation.v(t, Q, V)
    # and f_0(Q, V(Q)) = int.equation.f(t, Q, V, F)
    for i in 1:int.tableau.q.s
        tᵢ = int.params.t + int.Δt * int.tableau.q.c[i]
        int.equation.α(tᵢ, int.Q[i], int.ϑ[i])
        int.equation.v(tᵢ, int.Q[i], int.V[i])
        int.equation.g(tᵢ, int.Q[i], int.V[i], int.F[i])
    end

    # compute Jacobian of v via ForwardDiff
    for i in 1:int.tableau.q.s
        tᵢ = int.params.t + int.Δt * int.tableau.q.c[i]
        v_rev! = (v,q) -> int.equation.v(tᵢ,q,v)
        ForwardDiff.jacobian!(int.J[i], v_rev!, tV, int.Q[i])
    end

    # contract J with ϑ and add to G
    for l in 1:int.tableau.q.s
        for i in 1:int.equation.d
            int.G[l][i] = 0
            for j in 1:int.equation.d
                int.G[l][i] += int.ϑ[l][j] * int.J[l][j,i]
            end
        end
    end

    # solve linear system AP=δP for P

    # compute δP
    for l in 1:int.tableau.q.s
        for i in 1:int.equation.d
            # set δP = p
            δP[(l-1)*int.equation.d+i] = int.params.p[i]
            # add A(F+G) to δP
            for k in 1:int.tableau.q.s
                δP[(l-1)*int.equation.d+i] += int.Δt * int.tableau.q.a[l,k] * (int.F[k][i] + int.G[k][i])
            end
        end
    end

    # construct A = identity(sd×sd) + A ⊗ J
    for k in 1:int.tableau.q.s
        for l in 1:int.tableau.q.s
            for i in 1:int.equation.d
                for j in 1:int.equation.d
                    int.A[(k-1)*int.equation.d+i, (l-1)*int.equation.d+j] = int.Δt * int.tableau.q.a[k,l] * int.J[l][j,i]
                end
            end
        end
    end
    for i in 1:int.equation.d*int.tableau.q.s
        int.A[i,i] += one(DT)
    end

    # solve AP = δP and copy result to int.P
    lu = LUSolver(int.A, δP)
    factorize!(lu)
    solve!(lu)
    tP = reshape(lu.x, (int.equation.d, int.tableau.q.s))
    for l in 1:int.tableau.q.s
        int.P[l] .= tP[:,l]
    end

    # contract J with P and subtract from G, so that G = (ϑ-P)J
    for l in 1:int.tableau.q.s
        for i in 1:int.equation.d
            int.G[l][i] = 0
            for j in 1:int.equation.d
                int.G[l][i] += (int.ϑ[l][j] - int.P[l][j]) * int.J[l][j,i]
            end
        end
    end

    # println(int.ϑ)
    # println(int.P)
    # println(int.ϑ .- int.P)
    # println()

    # compute final update for p
    update_solution!(int.p[m], int.F, int.tableau.q.b, int.tableau.q.b̂, int.Δt)
    update_solution!(int.p[m], int.G, int.tableau.q.b, int.tableau.q.b̂, int.Δt)

    # copy solution to initial guess
    update!(int.iguess, m, sol.t[0] + n*int.Δt, int.q[m])

    # take care of periodic solutions
    cut_periodic_solution!(int.q[m], int.equation.periodicity)

    # copy to solution
    copy_solution!(sol, int.q[m], int.p[m], n, m)
end
