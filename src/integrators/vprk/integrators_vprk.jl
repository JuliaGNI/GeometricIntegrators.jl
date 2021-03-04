
"Parameters for right-hand side function of Variational Partitioned Runge-Kutta methods."
const ParametersVPRK = AbstractParametersVPRK{:vprk}

function IntegratorCache(params::ParametersVPRK{DT,TT,D,S}; kwargs...) where {DT,TT,D,S}
    IntegratorCacheVPRK{DT,D,S}(false; kwargs...)
end

function IntegratorCache{ST}(params::ParametersVPRK{DT,TT,D,S}; kwargs...) where {ST,DT,TT,D,S}
    IntegratorCacheVPRK{ST,D,S}(false; kwargs...)
end


"Variational partitioned Runge-Kutta integrator."
struct IntegratorVPRK{DT, TT, D, S, PT <: ParametersVPRK{DT,TT},
                                    ST <: NonlinearSolver{DT},
                                    IT <: InitialGuessIODE{TT}} <: AbstractIntegratorVPRK{DT,TT,D,S}
    params::PT
    solver::ST
    iguess::IT
    caches::CacheDict{PT}

    function IntegratorVPRK(params::ParametersVPRK{DT,TT,D,S}, solver::ST, iguess::IT, caches) where {DT,TT,D,S,ST,IT}
        new{DT, TT, D, S, typeof(params), ST, IT}(params, solver, iguess, caches)
    end

    function IntegratorVPRK{DT,D}(equations::NamedTuple, tableau::TableauVPRK{TT}, Δt::TT) where {DT,TT,D}
        # get number of stages
        S = tableau.s

        # create params
        params = ParametersVPRK{DT,D}(equations, tableau, Δt)

        # create cache dict
        caches = CacheDict(params)

        # create nonlinear solver
        solver = create_nonlinear_solver(DT, D*S, params, caches)

        # create initial guess
        iguess = InitialGuessIODE(get_config(:ig_interpolation), equations[:v̄], equations[:f̄], Δt)

        # create integrator
        IntegratorVPRK(params, solver, iguess, caches)
    end

    function IntegratorVPRK{DT,D}(ϑ::Function, f::Function, g::Function, v::Function,
                    tableau::TableauVPRK{TT}, Δt::TT; kwargs...) where {DT,TT,D}
        IntegratorVPRK{DT,D}(NamedTuple{(:ϑ,:f,:g,:v)}((ϑ, f, g, v)), tableau, Δt; kwargs...)
    end

    function IntegratorVPRK{DT,D}(ϑ::Function, f::Function, g::Function, v::Function, h::Function,
                    tableau::TableauVPRK{TT}, Δt::TT; kwargs...) where {DT,TT,D}
        IntegratorVPRK{DT,D}(NamedTuple{(:ϑ,:f,:g,:v,:h)}((ϑ, f, g, v, h)), tableau, Δt; kwargs...)
    end

    function IntegratorVPRK(equation::Union{IODE{DT},LODE{DT}}, tableau, Δt; kwargs...) where {DT}
        IntegratorVPRK{DT, ndims(equation)}(get_functions(equation), tableau, Δt; kwargs...)
    end
end


IntegratorVPRKpNone = IntegratorVPRK

function Integrators.get_internal_variables(int::IntegratorVPRK{DT,TT,D,S}) where {DT, TT, D, S}
    Q = create_internal_stage_vector(DT, D, S)
    P = create_internal_stage_vector(DT, D, S)
    V = create_internal_stage_vector(DT, D, S)
    F = create_internal_stage_vector(DT, D, S)

    solver = get_solver_status(int.solver)

    (Q=Q, P=P, V=V, F=F, solver=solver)
end


function initial_guess!(int::IntegratorVPRK{DT}, sol::AtomicSolutionPODE{DT},
                        cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT}
    for i in eachstage(int)
        evaluate!(int.iguess, sol.q̄, sol.p̄, sol.v̄, sol.f̄,
                              sol.q, sol.p, sol.v, sol.f,
                              cache.q̃, cache.ṽ,
                              tableau(int).q.c[i])

        for k in eachdim(int)
            int.solver.x[ndims(int)*(i-1)+k] = cache.ṽ[k]
        end
    end
end


"Compute stages of variational partitioned Runge-Kutta methods."
function Integrators.function_stages!(x::Vector{ST}, b::Vector{ST},
                params::ParametersVPRK{DT,TT,D,S},
                caches::CacheDict) where {ST,DT,TT,D,S}

    @assert length(x) == length(b)

    # get cache for internal stages
    cache = caches[ST]

    compute_stages!(x, cache.Q, cache.V, cache.P, cache.F, params)
    compute_rhs_vprk!(b, cache.P, cache.F, params)
    compute_rhs_vprk_correction!(b, cache.V, params)
end


function Integrators.integrate_step!(int::IntegratorVPRK{DT,TT}, sol::AtomicSolutionPODE{DT,TT},
                                     cache::IntegratorCacheVPRK{DT}=int.caches[DT]) where {DT,TT}
    # update nonlinear solver parameters from cache
    update_params!(int.params, sol)

    # compute initial guess
    initial_guess!(int, sol, cache)

    # reset solution
    reset!(sol, timestep(int))

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params)

    # compute vector fields at internal stages
    compute_stages!(int.solver.x, cache.Q, cache.V, cache.P, cache.F, int.params)

    # compute final update
    update_solution!(int, sol, cache)

    # copy solution to initial guess
    update_vector_fields!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f)

    # copy internal stage variables
    sol.internal[:Q] .= cache.Q
    sol.internal[:P] .= cache.P
    sol.internal[:V] .= cache.V
    sol.internal[:F] .= cache.F

    # copy solver status
    get_solver_status!(int.solver, sol.internal[:solver])
end
