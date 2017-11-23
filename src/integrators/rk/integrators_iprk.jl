
struct NonlinearFunctionCacheIPRK{ST}
    Q::Matrix{ST}
    V::Matrix{ST}
    P::Matrix{ST}
    F::Matrix{ST}
    Y::Matrix{ST}
    Z::Matrix{ST}

    v::Array{ST,1}
    f::Array{ST,1}
    y::Array{ST,1}
    z::Array{ST,1}

    function NonlinearFunctionCacheIPRK{ST}(D,S) where {ST}
        # create internal stage vectors
        Q = zeros(ST,D,S)
        V = zeros(ST,D,S)
        P = zeros(ST,D,S)
        F = zeros(ST,D,S)
        Y = zeros(ST,D,S)
        Z = zeros(ST,D,S)

        # create update vectors
        v = zeros(ST,D)
        f = zeros(ST,D)
        y = zeros(ST,D)
        z = zeros(ST,D)

        new(Q, V, P, F, Y, Z, v, f, y, z)
    end
end

"Parameters for right-hand side function of implicit partitioned Runge-Kutta methods."
mutable struct NonlinearFunctionParametersIPRK{DT,TT,VT,FT,D,S} <: NonlinearFunctionParameters{DT}
    f_v::VT
    f_f::FT

    Δt::TT

    t_q::CoefficientsRK{TT}
    t_p::CoefficientsRK{TT}

    t::TT

    q::Vector{DT}
    p::Vector{DT}

    function NonlinearFunctionParametersIPRK{DT,TT,VT,FT,D,S}(f_v, f_f, Δt, t_q, t_p) where {DT,TT,VT,FT,D,S}
        # create solution vectors
        q = zeros(DT,D)
        p = zeros(DT,D)

        new(f_v, f_f, Δt, t_q, t_p, 0, q, p)
    end
end

"Compute stages of implicit partitioned Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST}, params::NonlinearFunctionParametersIPRK{DT,TT,VT,FT,D,S}) where {ST,DT,TT,VT,FT,D,S}
    cache = NonlinearFunctionCacheIPRK{ST}(D, S)

    quote
        compute_stages!(x, $cache.Q, $cache.V, $cache.P, $cache.F, $cache.Y, $cache.Z, params)

        # compute b = - [(Y-AV), (Z-AF)]
        for i in 1:S
            for k in 1:D
                b[2*(D*(i-1)+k-1)+1] = - $cache.Y[k,i]
                b[2*(D*(i-1)+k-1)+2] = - $cache.Z[k,i]
                for j in 1:S
                    b[2*(D*(i-1)+k-1)+1] += params.t_q.a[i,j] * $cache.V[k,j]
                    b[2*(D*(i-1)+k-1)+2] += params.t_p.a[i,j] * $cache.F[k,j]
                end
            end
        end
    end
end


@generated function compute_stages!(x::Vector{ST}, Q::Matrix{ST}, V::Matrix{ST},
                                                   P::Matrix{ST}, F::Matrix{ST},
                                                   Y::Matrix{ST}, Z::Matrix{ST},
                                                   params::NonlinearFunctionParametersIPRK{DT,TT,VT,FT,D,S}) where {ST,DT,TT,VT,FT,D,S}
    tQ = zeros(ST,D)
    tV = zeros(ST,D)
    tP = zeros(ST,D)
    tF = zeros(ST,D)

    quote
        local tqᵢ::TT
        local tpᵢ::TT

        for i in 1:S
            for k in 1:D
                # copy y to Y and Z
                Y[k,i] = x[2*(D*(i-1)+k-1)+1]
                Z[k,i] = x[2*(D*(i-1)+k-1)+2]

                # compute Q and P
                Q[k,i] = params.q[k] + params.Δt * Y[k,i]
                P[k,i] = params.p[k] + params.Δt * Z[k,i]
            end

            # compute f(X)
            tqᵢ = params.t + params.Δt * params.t_q.c[i]
            tpᵢ = params.t + params.Δt * params.t_p.c[i]

            simd_copy_xy_first!($tQ, Q, i)
            simd_copy_xy_first!($tP, P, i)
            params.f_v(tqᵢ, $tQ, $tP, $tV)
            params.f_f(tpᵢ, $tQ, $tP, $tF)
            simd_copy_yx_first!($tV, V, i)
            simd_copy_yx_first!($tF, F, i)
        end
    end
end


"Implicit partitioned Runge-Kutta integrator."
struct IntegratorIPRK{DT, TT, VT, FT, SPT, ST, IT} <: Integrator{DT, TT}
    equation::PODE{DT,TT,VT,FT}
    tableau::TableauIPRK{TT}
    Δt::TT

    cache::NonlinearFunctionCacheIPRK{DT}
    params::SPT
    solver::ST
    iguess::InitialGuessPODE{DT,TT,VT,FT,IT}

    q::Vector{Vector{Double{DT}}}
    p::Vector{Vector{Double{DT}}}
end

function IntegratorIPRK(equation::PODE{DT,TT,VT,FT}, tableau::TableauIPRK{TT}, Δt::TT) where {DT,TT,VT,FT}
    D = equation.d
    M = equation.n
    S = tableau.s

    # create solution vector for internal stages / nonlinear solver
    x = zeros(DT, 2*D*S)

    # create solution vectors
    q = Array{Vector{Double{DT}}}(M)
    p = Array{Vector{Double{DT}}}(M)

    for i in 1:M
        q[i] = zeros(Double{DT},D)
        p[i] = zeros(Double{DT},D)
    end

    # create cache for internal stage vectors and update vectors
    cache = NonlinearFunctionCacheIPRK{DT}(D,S)

    # create params
    params = NonlinearFunctionParametersIPRK{DT,TT,VT,FT,D,S}(
                                                equation.v, equation.f,
                                                Δt, tableau.q, tableau.p)

    # create rhs function for nonlinear solver
    function_stages = (x,b) -> function_stages!(x, b, params)

    # create solver
    solver = get_config(:nls_solver)(x, function_stages)

    # create initial guess
    iguess = InitialGuessPODE(interpolation, equation, Δt)

    # create integrator
    IntegratorIPRK{DT, TT, VT, FT, typeof(params), typeof(solver), typeof(iguess.int)}(
                                        equation, tableau, Δt, cache, params, solver, iguess, q, p)
end


function initialize!(int::IntegratorIPRK, sol::SolutionPODE, m::Int)
    @assert m ≥ 1
    @assert m ≤ sol.ni

    # copy initial conditions from solution
    get_initial_conditions!(sol, int.q[m], int.p[m], m)

    # initialise initial guess
    initialize!(int.iguess, m, sol.t[0], int.q[m], int.p[m])
end

function initial_guess!(int::IntegratorIPRK, m::Int)
    for i in 1:int.tableau.q.s
        evaluate!(int.iguess, m, int.cache.y, int.cache.z, int.cache.v, int.cache.f, int.tableau.q.c[i], int.tableau.p.c[i])
        for k in 1:int.equation.d
            int.cache.V[k,i] = int.cache.v[k]
            int.cache.F[k,i] = int.cache.f[k]
        end
    end
    for i in 1:int.tableau.q.s
        for k in 1:int.equation.d
            int.solver.x[2*(int.equation.d*(i-1)+k-1)+1] = 0
            int.solver.x[2*(int.equation.d*(i-1)+k-1)+2] = 0
            for j in 1:int.tableau.q.s
                int.solver.x[2*(int.equation.d*(i-1)+k-1)+1] += int.tableau.q.a[i,j] * int.cache.V[k,j]
                int.solver.x[2*(int.equation.d*(i-1)+k-1)+2] += int.tableau.p.a[i,j] * int.cache.F[k,j]
            end
        end
    end
end

"Integrate ODE with implicit partitioned Runge-Kutta integrator."
function integrate_step!(int::IntegratorIPRK{DT,TT,VT,FT}, sol::SolutionPODE{DT,TT,N}, m::Int, n::Int) where {DT,TT,VT,FT,N}
    # set time for nonlinear solver
    int.params.t  = sol.t[0] + (n-1)*int.Δt
    int.params.q .= int.q[m]
    int.params.p .= int.p[m]

    # compute initial guess
    initial_guess!(int, m)

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params, n)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params, n)

    # compute vector fields at internal stages
    compute_stages!(int.solver.x, int.cache.Q, int.cache.V, int.cache.P, int.cache.F, int.cache.Y, int.cache.Z, int.params)

    # compute final update
    update_solution!(int.q[m], int.cache.V, int.tableau.q.b, int.tableau.q.b̂, int.Δt)
    update_solution!(int.p[m], int.cache.F, int.tableau.p.b, int.tableau.p.b̂, int.Δt)

    # take care of periodic solutions
    cut_periodic_solution!(int.q[m], int.equation.periodicity)

    # copy to solution
    copy_solution!(sol, int.q[m], int.p[m], n, m)
end
