
"Parameters for right-hand side function of fully implicit Runge-Kutta methods."
mutable struct ParametersMidpointImplicit{DT, TT, D, S, ET <: NamedTuple} <: Parameters{DT,TT}
    equs::ET
    tab::Tableau{TT}
    Δt::TT

    t::TT
    q::Vector{DT}
    θ::Vector{DT}

    function ParametersMidpointImplicit{DT,D}(equs::ET, tab::Tableau{TT}, Δt::TT) where {DT, TT, D, ET <: NamedTuple}
        new{DT, TT, D, tab.s, ET}(equs, tab, Δt, zero(TT), zeros(DT,D), zeros(DT,D))
    end
end


@doc raw"""
Fully implicit Runge-Kutta integrator cache.

### Fields

* `q̃`: initial guess of solution
* `ṽ`: initial guess of vector field
* `s̃`: holds shift due to periodicity of solution
* `Q`: internal stages of solution
* `V`: internal stages of vector field
* `Θ`: implicit function of internal stages
* `F`: vector field of implicit function
"""
struct IntegratorCacheMidpointImplicit{DT,D,S} <: ODEIntegratorCache{DT,D}
    q::Vector{DT}
    v::Vector{DT}
    θ::Vector{DT}

    Q::Vector{DT}
    V::Vector{DT}
    F::Vector{DT}

    function IntegratorCacheMidpointImplicit{DT,D,S}() where {DT,D,S}
        new(zeros(DT,D), zeros(DT,D), zeros(DT,D),
            zeros(DT,D), zeros(DT,D), zeros(DT,D))
    end
end

function IntegratorCache{ST}(params::ParametersMidpointImplicit{DT,TT,D,S}; kwargs...) where {ST,DT,TT,D,S}
    IntegratorCacheMidpointImplicit{ST,D,S}(; kwargs...)
end

@inline CacheType(ST, params::ParametersMidpointImplicit{DT,TT,D,S}) where {DT,TT,D,S} = IntegratorCacheMidpointImplicit{ST,D,S}


"Fully implicit Runge-Kutta integrator."
struct IntegratorMidpointImplicit{DT, TT, D, S, PT <: ParametersMidpointImplicit{DT,TT},
                                    ST <: NonlinearSolver{DT},
                                    IT <: InitialGuessODE{DT,TT}} <: AbstractIntegratorRK{DT,TT}
    params::PT
    solver::ST
    iguess::IT
    caches::CacheDict{PT}

    function IntegratorMidpointImplicit(params::ParametersMidpointImplicit{DT,TT,D,S}, solver::ST, iguess::IT, caches) where {DT,TT,D,S,ST,IT}
        new{DT, TT, D, S, typeof(params), ST, IT}(params, solver, iguess, caches)
    end

    function IntegratorMidpointImplicit{DT,D}(equations::NamedTuple, tableau::Tableau{TT}, Δt::TT) where {DT,TT,D}
        # get number of stages
        S = tableau.s

        @assert S == 1

        # create params
        params = ParametersMidpointImplicit{DT,D}(equations, tableau, Δt)

        # create cache dict
        caches = CacheDict(params)

        # create solver
        solver = create_nonlinear_solver(DT, D*S, params, caches)

        # create initial guess
        iguess = InitialGuessODE{DT,D}(get_config(:ig_interpolation), equations[:v̄], Δt)

        # create integrator
        IntegratorMidpointImplicit(params, solver, iguess, caches)
    end

    # function IntegratorMidpointImplicit{DT,D}(v::Function, tableau::Tableau{TT}, Δt::TT; kwargs...) where {DT,TT,D}
    #     IntegratorMidpointImplicit{DT,D}(NamedTuple{(:v,)}((v,)), tableau, Δt; kwargs...)
    # end

    # function IntegratorMidpointImplicit{DT,D}(v::Function, h::Function, tableau::Tableau{TT}, Δt::TT; kwargs...) where {DT,TT,D}
    #     IntegratorMidpointImplicit{DT,D}(NamedTuple{(:v,:h)}((v,h)), tableau, Δt; kwargs...)
    # end

    function IntegratorMidpointImplicit(equation::IODE{DT,TT}, tableau::Tableau{TT}, Δt::TT; kwargs...) where {DT,TT}
        IntegratorMidpointImplicit{DT, ndims(equation)}(get_function_tuple(equation), tableau, Δt; kwargs...)
    end
end


@inline Base.ndims(::IntegratorMidpointImplicit{DT,TT,D,S}) where {DT,TT,D,S} = D


function initialize!(int::IntegratorMidpointImplicit, sol::AtomicSolutionODE)
    sol.t̄ = sol.t - timestep(int)

    equations(int)[:v](sol.t, sol.q, sol.v)

    initialize!(int.iguess, sol.t, sol.q, sol.v,
                            sol.t̄, sol.q̄, sol.v̄)
end


function update_params!(int::IntegratorMidpointImplicit, sol::AtomicSolutionPODE)
    # set time for nonlinear solver and copy previous solution
    int.params.t  = sol.t
    int.params.q .= sol.q
    equations(int)[:ϑ](sol.t, sol.q, sol.v, int.params.θ)
end


function initial_guess!(int::IntegratorMidpointImplicit{DT,TT}, sol::AtomicSolutionPODE{DT,TT},
                        cache::IntegratorCacheMidpointImplicit{DT}=int.caches[DT]) where {DT,TT}

    # compute initial guess for solution
    evaluate!(int.iguess, sol.q̄, sol.v̄, sol.q, sol.v, cache.q, cache.v, one(TT))
    int.solver.x .= cache.q 
end


function compute_stages!(x::Vector{ST}, Q::Vector{ST}, V::Vector{ST}, F::Vector{ST}, 
                                        q::Vector{ST}, v::Vector{ST}, θ::Vector{ST},
                                        params::ParametersMidpointImplicit{DT,TT,D,S}) where {ST,DT,TT,D,S}

    # copy x to V
    q .= x

    # compute Q and V
    Q .= (q .+ params.q) ./ 2
    V .= (q .- params.q) ./ params.Δt

    # compute F = f(Q,V)
    params.equs[:f](params.t + params.Δt / 2, Q, V, F)
 
    # compute Θ = ϑ(q̄)
    params.equs[:ϑ](params.t + params.Δt, q, v, θ)
end

"Compute stages of fully implicit Runge-Kutta methods."
function function_stages!(x::Vector{ST}, b::Vector{ST}, params::ParametersMidpointImplicit{DT,TT,D,S},
                          caches::CacheDict) where {ST,DT,TT,D,S}
    # get cache for internal stages
    cache = caches[ST]

    # compute stages from nonlinear solver solution x
    compute_stages!(x, cache.Q, cache.V, cache.F, cache.q, cache.v, cache.θ, params)

    # compute b
    b .= cache.θ .- params.θ .- params.Δt .* cache.F
end


function integrate_step!(int::IntegratorMidpointImplicit{DT,TT}, sol::AtomicSolutionPODE{DT,TT},
                         cache::IntegratorCacheMidpointImplicit{DT}=int.caches[DT]) where {DT,TT}

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
    compute_stages!(int.solver.x, cache.Q, cache.V, cache.F, cache.q, cache.v, cache.θ, int.params)

    # compute final update
    sol.q .= cache.q
    sol.p .= cache.θ

    # compute vector field for initial guess
    equations(int)[:v](sol.t, sol.q, sol.v)
    # update_vector_fields!(int.iguess, sol.t, sol.q, sol.v)
end
