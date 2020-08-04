
"Parameters for right-hand side function of Variational Partitioned Runge-Kutta methods."
const ParametersVPRKpMidpoint = AbstractParametersVPRK{:vprk_pmidpoint}


"Variational partitioned Runge-Kutta integrator."
struct IntegratorVPRKpMidpoint{DT, TT, D, S,
                PT <: ParametersVPRKpMidpoint{DT,TT},
                ST <: NonlinearSolver{DT},
                IT <: InitialGuessIODE{DT,TT}} <: AbstractIntegratorVPRKwProjection{DT,TT,D,S}

    params::PT
    solver::ST
    iguess::IT
    caches::CacheDict{PT}

    function IntegratorVPRKpMidpoint(params::ParametersVPRKpMidpoint{DT,TT,D,S}, solver::ST, iguess::IT, caches) where {DT,TT,D,S,ST,IT}
        new{DT, TT, D, S, typeof(params), ST, IT}(params, solver, iguess, caches)
    end

    function IntegratorVPRKpMidpoint{DT,D}(equations::NamedTuple, tableau::TableauVPRK{TT}, Δt::TT) where {DT,TT,D}
        # get number of stages
        S = tableau.s

        # create params
        R = TT[1, tableau.R∞]
        params = ParametersVPRKpMidpoint{DT,D}(equations, tableau, Δt, NamedTuple{(:R,)}((R,)))

        # create cache dict
        caches = CacheDict(params)

        # create solver
        solver = create_nonlinear_solver(DT, D*(S+2), params, caches)

        # create initial guess
        iguess = InitialGuessIODE{DT,D}(get_config(:ig_interpolation), equations[:v], equations[:f], Δt)

        # create integrator
        IntegratorVPRKpMidpoint(params, solver, iguess, caches)
    end

    function IntegratorVPRKpMidpoint(equation::IODE{DT,TT}, tableau::TableauVPRK{TT}, Δt::TT; kwargs...) where {DT,TT}
        IntegratorVPRKpMidpoint{DT, ndims(equation)}(get_function_tuple(equation), tableau, Δt; kwargs...)
    end
end


function initial_guess!(int::IntegratorVPRKpMidpoint{DT,TT}, sol::AtomicSolutionPODE{DT,TT},
                        cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT,TT}
    for i in eachstage(int)
        evaluate!(int.iguess, sol.q̅, sol.p̅, sol.v̅, sol.f̅,
                              sol.q, sol.p, sol.v, sol.f,
                              cache.q̃, cache.ṽ,
                              tableau(int).q.c[i])

        for k in eachdim(int)
            int.solver.x[ndims(int)*(i-1)+k] = cache.ṽ[k]
        end
    end

    evaluate!(int.iguess, sol.q, sol.p, sol.v, sol.f,
                          sol.q̅, sol.p̅, sol.v̅, sol.f̅,
                          cache.q̃, one(TT))
    for k in eachdim(int)
        int.solver.x[ndims(int)*(nstages(int)+0)+k] = cache.q̃[k]
    end
    for k in eachdim(int)
        int.solver.x[ndims(int)*(nstages(int)+1)+k] = 0
    end
end


function compute_projection_vprk!(x::Vector{ST},
                q::SolutionVector{ST}, p::SolutionVector{ST}, λ::SolutionVector{ST},
                V::Vector{Vector{ST}}, U::Vector{Vector{ST}}, G::Vector{Vector{ST}},
                params::ParametersVPRKpMidpoint{DT,TT,D,S}) where {ST,DT,TT,D,S}

    # create temporary variables
    local t₀::TT = params.t̅
    local t₁::TT = params.t̅ + params.Δt
    local tₘ::TT = (t₀+t₁)/2
    local y::ST
    local q̃ = zeros(ST,D)

    # copy x to λ and q̅
    for k in 1:D
        q[k] = x[D*(S+0)+k]
        λ[k] = x[D*(S+1)+k]
    end

    # compute U=λ
    U[1] .= λ
    U[2] .= λ

    # compute G=g(q,λ)
    for k in 1:D
        y = 0
        for j in 1:S
            y += params.tab.q.b[j] * V[j][k]
        end
        q̃[k] = params.q̅[k] + 0.5 * params.Δt * y + params.Δt * params.pparams[:R][1] * U[1][k]
        # qm[k] = params.q̅[k] + 0.5 * params.Δt * y + 0.5 * params.Δt * params.R[1] * U[k,1] + 0.5 * params.Δt * params.R[2] * U[k,2]
    end

    # println("q̃mid = ", q̃)
    # println("qmid = ", qm)

    params.equ[:g](tₘ, q̃, λ, G[1])
    G[2] .= G[1]

    # compute p=ϑ(q)
    params.equ[:ϑ](t₁, q, λ, p)
end


"Compute stages of variational partitioned Runge-Kutta methods."
function Integrators.function_stages!(x::Vector{ST}, b::Vector{ST},
                params::ParametersVPRKpMidpoint{DT,TT,D,S},
                caches::CacheDict) where {ST,DT,TT,D,S}

    # get cache for internal stages
    cache = caches[ST]

    compute_stages!(x, cache.q̃, cache.p̃, cache.λ, cache.Q, cache.V, cache.U, cache.P, cache.F, cache.G, params)

    # compute b = - [P-AF-U]
    compute_rhs_vprk!(b, cache.P, cache.F, cache.G, params)

    # compute b = - [p-bF-G]
    compute_rhs_vprk_projection_p!(b, cache.p̃, cache.F, cache.G, D*(S+0), params)

    # compute b = - [q-bV-U]
    compute_rhs_vprk_projection_q!(b, cache.q̃, cache.V, cache.U, D*(S+1), params)

    compute_rhs_vprk_correction!(b, cache.V, params)
end


"Integrate ODE with variational partitioned Runge-Kutta integrator."
function Integrators.integrate_step!(int::IntegratorVPRKpMidpoint{DT,TT}, sol::AtomicSolutionPODE{DT,TT},
                                     cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT,TT}
    # update nonlinear solver parameters from cache
    update_params!(int.params, sol)

    # compute initial guess
    initial_guess!(int, sol, cache)

    # reset cache
    reset!(sol, timestep(int))

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params)

    # compute vector fields at internal stages and projection vector fields
    compute_stages!(int.solver.x,
                    cache.q̃, cache.p̃, cache.λ,
                    cache.Q, cache.V, cache.U,
                    cache.P, cache.F, cache.G, int.params)

    # compute unprojected solution
    update_solution!(int, sol, cache)

    # add projection to solution
    project_solution!(int, sol, int.params.pparams[:R], cache)

    # copy solution to initial guess
    update_vector_fields!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f)
end
