"""
Implicit Euler Method.

$(reference(Val(:ImplicitEuler)))
"""
struct ImplicitEuler <: ODEMethod end


@doc raw"""
Implicit Euler integrator cache.
"""
struct IntegratorCacheImplicitEuler{DT,D} <: ODEIntegratorCache{DT,D}
    x::Vector{DT}
    q̄::Vector{DT}
    q::Vector{DT}
    v::Vector{DT}
    v̄::Vector{DT}

    function IntegratorCacheImplicitEuler{DT,D}() where {DT,D}
        x = zeros(DT, D)
        q̄ = zeros(DT, D)
        q = zeros(DT, D)
        v = zeros(DT, D)
        v̄ = zeros(DT, D)
        new(x, q̄, q, v, v̄)
    end
end

nlsolution(cache::IntegratorCacheImplicitEuler) = cache.x
reset!(cache::IntegratorCacheImplicitEuler, t, q, λ = missing) = copyto!(cache.q̄, q)

function Cache{ST}(problem::AbstractProblem, method::ImplicitEuler; kwargs...) where {ST}
    IntegratorCacheImplicitEuler{ST, ndims(problem)}(; kwargs...)
end

@inline CacheType(ST, problem::AbstractProblem, method::ImplicitEuler) = IntegratorCacheImplicitEuler{ST, ndims(problem)}


const IntegratorImplicitEuler{DT,TT} = Integrator{<:Union{ODEProblem{DT,TT}, DAEProblem{DT,TT}, SubstepProblem{DT,TT}}, <:ImplicitEuler}

solversize(problem::Union{ODEProblem, DAEProblem, SubstepProblem}, ::ImplicitEuler) = ndims(problem)

default_solver(::ImplicitEuler) = Newton()
default_iguess(::ImplicitEuler) = HermiteExtrapolation()


function initial_guess!(
    solstep::Union{SolutionStepODE{DT}, SolutionStepDAE{DT}, SubstepProblem{DT}},
    problem::Union{ODEProblem, DAEProblem, SubstepProblem},
    ::ImplicitEuler, 
    caches::CacheDict, 
    ::NonlinearSolver, 
    iguess::Union{InitialGuess,Extrapolation}) where {DT}
    
    # obtain cache
    cache = caches[DT]

    # compute initial guess
    initialguess!(solstep.t, cache.q, cache.v, solstep, problem, iguess)

    # assemble initial guess for nonlinear solver solution vector
    nlsolution(cache) .= cache.v
end

function components!(x::AbstractVector{ST}, int::IntegratorImplicitEuler) where {ST}
    local q = cache(int, ST).q
    local q̄ = cache(int, ST).q̄
    local v = cache(int, ST).v
    local v̄ = cache(int, ST).v̄

    # compute q = q̄ + Δt * x (v = x)
    v̄ .= x
    q .= q̄ .+ timestep(int) .* v̄

    # compute v = v(q)
    equations(int).v(v, solstep(int).t, q)
end


function residual!(b::AbstractVector{ST}, int::IntegratorImplicitEuler) where {ST}
    # temporary variables
    local y::ST

    # get cache for internal stages
    local v = cache(int, ST).v
    local v̄ = cache(int, ST).v̄

    # compute b = - (v-v)
    b .= v .- v̄
end


# Compute stages of implicit Euler methods.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, int::IntegratorImplicitEuler) where {ST}
    @assert axes(x) == axes(b)

    # copy previous solution from solstep to cache
    reset!(cache(int, ST), current(solstep(int))...)

    # compute stages from nonlinear solver solution x
    components!(x, int)

    # compute residual vector
    residual!(b, int)
end


function update!(x::AbstractVector{DT}, int::IntegratorImplicitEuler) where {DT}
    # copy previous solution from solstep to cache
    reset!(cache(int, DT), current(solstep(int))...)

    # compute vector field at internal stages
    components!(x, int)

    # compute final update
    update!(solstep(int), cache(int, DT).v, timestep(int))
end


function integrate_step!(int::IntegratorImplicitEuler)
    # call nonlinear solver
    solve!(nlsolution(int), (b,x) -> residual!(b, x, int), solver(int))

    # print solver status
    # println(status(solver))

    # check if solution contains NaNs or error bounds are violated
    # println(meets_stopping_criteria(status(solver)))

    # compute final update
    update!(nlsolution(int), int)

    # update vector field for initial guess
    update_vector_fields!(solstep(int), problem(int))
end
