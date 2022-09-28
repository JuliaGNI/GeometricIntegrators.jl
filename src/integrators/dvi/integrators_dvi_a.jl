
import ForwardDiff

"Parameters for right-hand side function of degenerate variational integrator."
mutable struct ParametersDVIA{DT, TT, D, ET <: NamedTuple} <: Parameters{DT,TT}
    equs::ET
    Δt::TT

    t::TT
    q::Vector{DT}
    v::Vector{DT}
    θ::Vector{DT}

    function ParametersDVIA{DT,D}(equs::ET, Δt::TT) where {DT, TT, D, ET <: NamedTuple}
        new{DT, TT, D, ET}(equs, Δt, zero(TT), zeros(DT,D), zeros(DT,D), zeros(DT,D))
    end
end


@doc raw"""
Degenerate variational integrator cache.

### Fields

* `q`: internal stages of solution
* `v`: internal stages of vector field
* `Θ`: implicit function evaluated on solution
* `f`: vector field of implicit function
"""
struct IntegratorCacheDVIA{DT,D,S} <: ODEIntegratorCache{DT,D}
    q::Vector{DT}
    v::Vector{DT}
    θ::Vector{DT}
    f::Vector{DT}

    v̄::Vector{DT}
    θ̄::Vector{DT}
    f̄::Vector{DT}

    function IntegratorCacheDVIA{DT,D,S}() where {DT,D,S}
        new(zeros(DT,D), zeros(DT,D), zeros(DT,D), zeros(DT,D),
                         zeros(DT,D), zeros(DT,D), zeros(DT,D))
    end
end

function IntegratorCache{ST}(params::ParametersDVIA{DT,TT,D,S}; kwargs...) where {ST,DT,TT,D,S}
    IntegratorCacheDVIA{ST,D,S}(; kwargs...)
end

@inline CacheType(ST, ::ParametersDVIA{DT,TT,D,S}) where {DT,TT,D,S} = IntegratorCacheDVIA{ST,D,S}


"""
Implicit midpoint Runge-Kutta integrator.
"""
struct IntegratorDVIA{DT, TT, D, PT <: ParametersDVIA{DT,TT},
                                 ST <: NonlinearSolver{DT},
                                 IT <: InitialGuessIODE{TT}} <: AbstractIntegratorRK{DT,TT}
    params::PT
    solver::ST
    iguess::IT
    caches::CacheDict{PT}

    function IntegratorDVIA(params::ParametersDVIA{DT,TT,D}, solver::ST, iguess::IT, caches) where {DT,TT,D,ST,IT}
        new{DT, TT, D, typeof(params), ST, IT}(params, solver, iguess, caches)
    end

    function IntegratorDVIA{DT,D}(equations::NamedTuple, Δt::TT) where {DT,TT,D}
        # create params
        params = ParametersDVIA{DT,D}(equations, Δt)

        # create cache dict
        caches = CacheDict(params)

        # create solver
        solver = create_nonlinear_solver(DT, 2D, params, caches)

        # create initial guess
        iguess = InitialGuessIODE(get_config(:ig_extrapolation), equations[:v̄], equations[:f̄], Δt)

        # create integrator
        IntegratorDVIA(params, solver, iguess, caches)
    end

    function IntegratorDVIA(problem::Union{IODEProblem{DT}, LODEProblem{DT}}; kwargs...) where {DT,TT}
        IntegratorDVIA{DT, ndims(problem)}(functions(problem), timestep(problem); kwargs...)
    end
end


@inline Base.ndims(::IntegratorDVIA{DT,TT,D}) where {DT,TT,D} = D


function initialize!(int::IntegratorDVIA, sol::SolutionStepPODE)
    sol.t̄ = sol.t - timestep(int)

    equation(int, :v̄)(sol.t, sol.q, sol.v)
    equation(int, :f̄)(sol.t, sol.q, sol.v, sol.f)
    equation(int, :ϑ)(sol.t, sol.q, sol.v, sol.p)

    initialize!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f,
                            sol.t̄, sol.q̄, sol.p̄, sol.v̄, sol.f̄)
end


function update_params!(int::IntegratorDVIA, sol::SolutionStepPODE)
    # set time for nonlinear solver and copy previous solution
    int.params.t  = sol.t
    int.params.q .= sol.q
    int.params.v .= sol.v
    equations(int)[:ϑ](sol.t, sol.q, sol.v, int.params.θ)
end


function initial_guess!(int::IntegratorDVIA{DT,TT,D}, sol::SolutionStepPODE{DT,TT},
                        cache::IntegratorCacheDVIA{DT}=int.caches[DT]) where {DT,TT,D}

    # compute initial guess for solution
    evaluate!(int.iguess, sol.q̄, sol.p̄, sol.v̄, sol.f̄,
                          sol.q, sol.p, sol.v, sol.f,
                          cache.q, cache.v, one(TT))

    int.solver.x[1:D] .= cache.q

    for k in 1:div(D,2)
        int.solver.x[D+k] = cache.v[k]
        int.solver.x[D+div(D,2)+k] = sol.v[k]
    end
end


function compute_stages!(x::Vector{ST}, q::Vector{ST}, v::Vector{ST}, θ::Vector{ST}, f::Vector{ST},
                                                       v̄::Vector{ST}, θ̄::Vector{ST}, f̄::Vector{ST},
                                                       params::ParametersDVIA{DT,TT,D}) where {ST,DT,TT,D}

    # set some local variables for convenience and clarity
    local t = params.t + params.Δt
    local t̄ = params.t
    local q̄ = params.q
    θ̄  .= params.θ
    
    # copy x to q
    q .= x[1:D]

    # copy x to v and v̄
    for k in 1:div(D,2)
        v[k] = x[D+k]
        v̄[k] = x[D+div(D,2)+k]
        v[div(D,2)+k] = 0
        v̄[div(D,2)+k] = 0
    end

    # compute f = f(q,v)
    params.equs[:f](t, q, v, f)
    params.equs[:f](t̄, q̄, v̄, f̄)
 
    # compute Θ = ϑ(q,v)
    params.equs[:ϑ](t, q, v, θ)
    # params.equs[:ϑ](t̄, q̄, v̄, θ̄)
end


function function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersDVIA{DT,TT,D},
                          caches::CacheDict) where {ST,DT,TT,D}
    # get cache for internal stages
    cache = caches[ST]

    # compute stages from nonlinear solver solution x
    compute_stages!(x, cache.q, cache.v, cache.θ, cache.f, cache.v̄, cache.θ̄, cache.f̄, params)

    # compute b
    for k in 1:div(D,2)
        b[k]   = cache.θ[k] - cache.θ̄[k]  - params.Δt * cache.f̄[k]
        b[D+k] = cache.q[k] - params.q[k] - params.Δt * cache.v[k]
    end

    for k in div(D,2)+1:D
        b[k]   = params.Δt * cache.f̄[k]
        b[D+k] = params.Δt * cache.f[k]
    end
end


function integrate_step!(int::IntegratorDVIA{DT,TT}, sol::SolutionStepPODE{DT,TT},
                         cache::IntegratorCacheDVIA{DT}=int.caches[DT]) where {DT,TT}

    # update nonlinear solver parameters from atomic solution
    update_params!(int, sol)

    # compute initial guess
    initial_guess!(int, sol, cache)

    # reset atomic solution
    reset!(sol)

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params)

    # compute vector field at internal stages
    compute_stages!(int.solver.x, cache.q, cache.v, cache.θ, cache.f, cache.v̄, cache.θ̄, cache.f̄, int.params)

    # compute final update
    sol.q .= cache.q
    sol.v .= cache.v
    sol.p .= cache.θ
    sol.f .= cache.f
end
