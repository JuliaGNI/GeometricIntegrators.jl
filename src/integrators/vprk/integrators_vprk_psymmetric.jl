
"Parameters for right-hand side function of variational partitioned Runge-Kutta methods."
mutable struct ParametersVPRKpSymmetric{DT, TT, ET <: IODE{DT,TT}, D, S} <: AbstractParametersVPRK{DT,TT,ET,D,S}
    equ::ET
    tab::TableauVPRK{TT}
    Δt::TT

    R::Vector{TT}

    t̅::TT
    q̅::Vector{DT}
    p̅::Vector{DT}
end

function ParametersVPRKpSymmetric(equ::ET, tab::TableauVPRK{TT}, Δt::TT) where {DT, TT, ET <: IODE{DT,TT}}
    R = convert(Vector{TT}, [1, tab.R∞])
    q̅ = zeros(DT, ndims(equ))
    p̅ = zeros(DT, ndims(equ))

    ParametersVPRKpSymmetric{DT, TT, ET, ndims(equ), tab.s}(equ, tab, Δt, R, zero(TT), q̅, p̅)
end


"Variational partitioned Runge-Kutta integrator."
struct IntegratorVPRKpSymmetric{DT, TT, PT <: ParametersVPRKpSymmetric{DT,TT},
                                        ST <: NonlinearSolver{DT},
                                        IT <: InitialGuessPODE{DT,TT}, D, S} <: AbstractIntegratorVPRKwProjection{DT,TT}

    params::PT
    solver::ST
    iguess::IT
    cache::IntegratorCacheVPRK{DT,D,S}
end

function IntegratorVPRKpSymmetric(equation::ET, tableau::TableauVPRK{TT}, Δt::TT) where {DT, TT, ET <: IODE{DT,TT}}
    D = equation.d
    M = equation.n
    S = tableau.s

    # create params
    params = ParametersVPRKpSymmetric(equation, tableau, Δt)

    # create solver
    solver = create_nonlinear_solver(DT, D*(S+2), params)

    # create initial guess
    iguess = InitialGuessPODE(get_config(:ig_interpolation), equation, Δt)

    # create cache
    cache = IntegratorCacheVPRK{DT,D,S}(true)

    # create integrator
    IntegratorVPRKpSymmetric{DT, TT, typeof(params), typeof(solver), typeof(iguess), D, S}(params, solver, iguess, cache)
end


@inline Base.ndims(int::IntegratorVPRKpSymmetric{DT,TT,PT,ST,IT,D,S}) where {DT,TT,PT,ST,IT,D,S} = D


function compute_projection_vprk!(x::Vector{ST},
                q::SolutionVector{ST}, p::SolutionVector{ST}, λ::SolutionVector{ST},
                V::Vector{Vector{ST}}, U::Vector{Vector{ST}}, G::Vector{Vector{ST}},
                params::ParametersVPRKpSymmetric{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

    local t₀::TT = params.t̅
    local t₁::TT = params.t̅ + params.Δt

    # copy x to λ and q̅
    for k in 1:D
        q[k] = x[D*(S+0)+k]
        λ[k] = x[D*(S+1)+k]
    end

    # compute U=λ
    U[1] .= λ
    U[2] .= λ

    # compute G=g(q,λ)
    params.equ.g(t₀, params.q̅, λ, G[1])
    params.equ.g(t₁, q, λ, G[2])

    # compute p=ϑ(q)⋅p
    params.equ.ϑ(t₁, q, λ, p)
end


"Compute stages of variational partitioned Runge-Kutta methods."
@generated function function_stages!(x::Vector{ST}, b::Vector{ST},
                params::ParametersVPRKpSymmetric{DT,TT,ET,D,S}
            ) where {ST,DT,TT,ET,D,S}

    cache = IntegratorCacheVPRK{ST, D, S}(true)

    quote
        compute_stages!(x, $cache.q̃, $cache.p̃, $cache.λ, $cache.Q, $cache.V, $cache.U, $cache.P, $cache.F, $cache.G, params)

        # compute b = - [P-AF-U]
        compute_rhs_vprk!(b, $cache.P, $cache.F, $cache.G, params)

        # compute b = - [p-bF-G]
        compute_rhs_vprk_projection_p!(b, $cache.p̃, $cache.F, $cache.G, D*(S+0), params)

        # compute b = - [q-bV-U]
        compute_rhs_vprk_projection_q!(b, $cache.q̃, $cache.V, $cache.U, D*(S+1), params)

        compute_rhs_vprk_correction!(b, $cache.V, params)
    end
end


function initial_guess!(int::IntegratorVPRKpSymmetric{DT,TT}, sol::AtomicSolutionPODE{DT,TT}) where {DT,TT}
    for i in eachstage(int)
        evaluate!(int.iguess, sol.q, sol.p, sol.v, sol.f,
                              sol.q̅, sol.p̅, sol.v̅, sol.f̅,
                              int.cache.q̃, int.cache.ṽ,
                              tableau(int).q.c[i])

        for k in eachdim(int)
            int.solver.x[ndims(int)*(i-1)+k] = int.cache.ṽ[k]
        end
    end

    evaluate!(int.iguess, sol.q, sol.p, sol.v, sol.f,
                          sol.q̅, sol.p̅, sol.v̅, sol.f̅,
                          int.cache.q̃, one(TT))
    for k in eachdim(int)
        int.solver.x[ndims(int)*(nstages(int)+0)+k] = int.cache.q̃[k]
    end
    for k in eachdim(int)
        int.solver.x[ndims(int)*(nstages(int)+1)+k] = 0
    end
end

"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate_step!(int::IntegratorVPRKpSymmetric{DT,TT}, sol::AtomicSolutionPODE{DT,TT}) where {DT,TT}
    # update nonlinear solver parameters from cache
    update_params!(int.params, sol)

    # compute initial guess
    initial_guess!(int, sol)

    # reset solution
    reset!(sol, timestep(int))

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params)

    # compute vector fields at internal stages and projection vector fields
    compute_stages!(int.solver.x,
                    int.cache.q̃, int.cache.p̃, int.cache.λ,
                    int.cache.Q, int.cache.V, int.cache.U,
                    int.cache.P, int.cache.F, int.cache.G, int.params)

    # compute unprojected solution
    update_solution!(int, sol)

    # add projection to solution
    project_solution!(int, sol, int.params.R)

    # copy solution to initial guess
    update!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f)
end
