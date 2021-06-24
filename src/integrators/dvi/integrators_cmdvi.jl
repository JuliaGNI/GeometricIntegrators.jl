
import ForwardDiff

"Parameters for right-hand side function of degenerate variational integrator."
mutable struct ParametersCMDVI{DT, TT, D, ET <: NamedTuple} <: Parameters{DT,TT}
    equs::ET
    Δt::TT

    t::TT
    q::Vector{DT}
    v::Vector{DT}
    θ::Vector{DT}

    function ParametersCMDVI{DT,D}(equs::ET, Δt::TT) where {DT, TT, D, ET <: NamedTuple}
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
struct IntegratorCacheCMDVI{DT,D,S} <: ODEIntegratorCache{DT,D}
    q::Vector{DT}
    v::Vector{DT}
    θ::Vector{DT}

    q̄::Vector{DT}
    θ̄::Vector{DT}
    f̄::Vector{DT}

    function IntegratorCacheCMDVI{DT,D,S}() where {DT,D,S}
        new(zeros(DT,D), zeros(DT,D), zeros(DT,D),
            zeros(DT,D), zeros(DT,D), zeros(DT,D))
    end
end

function IntegratorCache{ST}(params::ParametersCMDVI{DT,TT,D,S}; kwargs...) where {ST,DT,TT,D,S}
    IntegratorCacheCMDVI{ST,D,S}(; kwargs...)
end

@inline CacheType(ST, ::ParametersCMDVI{DT,TT,D,S}) where {DT,TT,D,S} = IntegratorCacheCMDVI{ST,D,S}


"""
Implicit midpoint Runge-Kutta integrator.
"""
struct IntegratorCMDVI{DT, TT, D, PT <: ParametersCMDVI{DT,TT},
                                   ST <: NonlinearSolver{DT},
                                   IT <: InitialGuessIODE{TT}} <: AbstractIntegratorRK{DT,TT}
    params::PT
    solver::ST
    iguess::IT
    caches::CacheDict{PT}

    function IntegratorCMDVI(params::ParametersCMDVI{DT,TT,D}, solver::ST, iguess::IT, caches) where {DT,TT,D,ST,IT}
        new{DT, TT, D, typeof(params), ST, IT}(params, solver, iguess, caches)
    end

    function IntegratorCMDVI{DT,D}(equations::NamedTuple, Δt::TT) where {DT,TT,D}
        # create params
        params = ParametersCMDVI{DT,D}(equations, Δt)

        # create cache dict
        caches = CacheDict(params)

        # create solver
        solver = create_nonlinear_solver(DT, 2D, params, caches)

        # create initial guess
        iguess = InitialGuessIODE(get_config(:ig_extrapolation), equations[:v̄], equations[:f̄], Δt)

        # create integrator
        IntegratorCMDVI(params, solver, iguess, caches)
    end

    function IntegratorCMDVI(equation::Union{IODE{DT}, LODE{DT}}, Δt::TT; kwargs...) where {DT,TT}
        IntegratorCMDVI{DT, ndims(equation)}(get_functions(equation), Δt; kwargs...)
    end
end


@inline Base.ndims(::IntegratorCMDVI{DT,TT,D}) where {DT,TT,D} = D


function initialize!(int::IntegratorCMDVI, sol::AtomicSolutionPODE)
    sol.t̄ = sol.t - timestep(int)

    equation(int, :v̄)(sol.t, sol.q, sol.v)
    equation(int, :f̄)(sol.t, sol.q, sol.v, sol.f)
    equation(int, :ϑ)(sol.t, sol.q, sol.v, sol.p)

    initialize!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f,
                            sol.t̄, sol.q̄, sol.p̄, sol.v̄, sol.f̄)
end


function update_params!(int::IntegratorCMDVI, sol::AtomicSolutionPODE)
    # set time for nonlinear solver and copy previous solution
    int.params.t  = sol.t
    int.params.q .= sol.q
    int.params.v .= sol.v
    equations(int)[:ϑ](sol.t, sol.q, sol.v, int.params.θ)
end


function initial_guess!(int::IntegratorCMDVI{DT,TT,D}, sol::AtomicSolutionPODE{DT,TT},
                        cache::IntegratorCacheCMDVI{DT}=int.caches[DT]) where {DT,TT,D}

    # compute initial guess for solution
    evaluate!(int.iguess, sol.q̄, sol.p̄, sol.v̄, sol.f̄,
                          sol.q, sol.p, sol.v, sol.f,
                          cache.q, cache.v, one(TT))

    int.solver.x[1:D] .= cache.q                           # q(n+1)

    # compute initial guess for solution
    evaluate!(int.iguess, sol.q̄, sol.p̄, sol.v̄, sol.f̄,
                          sol.q, sol.p, sol.v, sol.f,
                          cache.q, cache.v, one(TT)/2)

    offset_v = D
    offset_x = D + div(D,2)
    for k in 1:div(D,2)
        int.solver.x[offset_v+k] = cache.v[k]              # v¹(n+1/2)
        int.solver.x[offset_x+k] = cache.q[div(D,2)+k]     # q²(n+1/2)
    end
end


function compute_stages!(x::Vector{ST}, q::Vector{ST}, v::Vector{ST}, θ::Vector{ST},
                                        q̄::Vector{ST}, θ̄::Vector{ST}, f̄::Vector{ST},
                                        params::ParametersCMDVI{DT,TT,D}) where {ST,DT,TT,D}

    # set some local variables for convenience and clarity
    local t⁻ = params.t
    local t̄  = params.t + params.Δt / 2
    local t⁺ = params.t + params.Δt
    
    # copy x to q
    q .= x[1:D]

    # copy x to q⁻, q⁺ and v
    for k in 1:div(D,2)
        q̄[k] = (params.q[k] + q[k]) / 2
        q̄[div(D,2)+k] = x[D+div(D,2)+k]

        v[k]          = x[D+k]
        v[div(D,2)+k] = 0
    end

    # compute f = f(q,v)
    params.equs[:f](t̄, q̄, v, f̄)
 
    # compute Θ = ϑ(q,v)
    params.equs[:ϑ](t̄, q̄, v, θ̄)
    params.equs[:ϑ](t⁺, q, v, θ)
end


function function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersCMDVI{DT,TT,D},
                          caches::CacheDict) where {ST,DT,TT,D}
    # get cache for internal stages
    cache = caches[ST]

    # compute stages from nonlinear solver solution x
    compute_stages!(x, cache.q, cache.v, cache.θ, cache.q̄, cache.θ̄, cache.f̄, params)

    # compute b
    b[1:D] .= cache.θ̄  .- params.θ .- params.Δt .* cache.f̄ ./ 2
    
    for k in 1:div(D,2)
        b[D+k] = cache.q[k] - params.q[k] - params.Δt * cache.v[k]
        b[D+div(D,2)+k] = cache.θ[k] - params.θ[k] - params.Δt * cache.f̄[k]
    end
end


function integrate_step!(int::IntegratorCMDVI{DT,TT}, sol::AtomicSolutionPODE{DT,TT},
                         cache::IntegratorCacheCMDVI{DT}=int.caches[DT]) where {DT,TT}

    # update nonlinear solver parameters from atomic solution
    update_params!(int, sol)

    # compute initial guess
    initial_guess!(int, sol, cache)

    # reset atomic solution
    reset!(sol, timestep(int))

    # call nonlinear solver
    solve!(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params)

    # compute vector field at internal stages
    compute_stages!(int.solver.x, cache.q, cache.v, cache.θ, cache.q̄, cache.θ̄, cache.f̄, int.params)

    # compute final update
    sol.q .= cache.q
    sol.p .= cache.θ

    # copy solution to initial guess
    update_vector_fields!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f)
end
