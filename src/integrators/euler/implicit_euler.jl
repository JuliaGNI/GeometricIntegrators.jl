"""
Implicit Euler Method.

$(reference(Val(:ImplicitEuler)))
"""
struct ImplicitEuler <: ODEMethod end

isexplicit(method::ImplicitEuler) = false
isimplicit(method::ImplicitEuler) = true
issymmetric(method::ImplicitEuler) = false
issymplectic(method::ImplicitEuler) = false


@doc raw"""
Implicit Euler integrator cache.
"""
struct ImplicitEulerCache{DT,D} <: ODEIntegratorCache{DT,D}
    x::Vector{DT}
    q̄::Vector{DT}
    q::Vector{DT}
    v::Vector{DT}
    v̄::Vector{DT}

    function ImplicitEulerCache{DT,D}() where {DT,D}
        x = zeros(DT, D)
        q̄ = zeros(DT, D)
        q = zeros(DT, D)
        v = zeros(DT, D)
        v̄ = zeros(DT, D)
        new(x, q̄, q, v, v̄)
    end
end

nlsolution(cache::ImplicitEulerCache) = cache.x

reset!(cache::ImplicitEulerCache, t, q, λ = missing) = copyto!(cache.q̄, q)

function Cache{ST}(problem::AbstractProblem, method::ImplicitEuler; kwargs...) where {ST}
    ImplicitEulerCache{ST, ndims(problem)}(; kwargs...)
end

@inline CacheType(ST, problem::AbstractProblem, method::ImplicitEuler) = ImplicitEulerCache{ST, ndims(problem)}


solversize(problem::AbstractProblemODE, ::ImplicitEuler) = ndims(problem)

default_solver(::ImplicitEuler) = Newton()
default_iguess(::ImplicitEuler) = HermiteExtrapolation()


function initial_guess!(int::GeometricIntegrator{<:ImplicitEuler})
    # compute initial guess
    initialguess!(solstep(int).t, cache(int).q, cache(int).v, solstep(int), problem(int), iguess(int))

    # assemble initial guess for nonlinear solver solution vector
    nlsolution(int) .= cache(int).v
end

function components!(x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:ImplicitEuler}) where {ST}
    local q = cache(int, ST).q
    local q̄ = cache(int, ST).q̄
    local v = cache(int, ST).v
    local v̄ = cache(int, ST).v̄

    # compute q = q̄ + Δt * x (v = x)
    v̄ .= x
    q .= q̄ .+ timestep(int) .* v̄

    # compute v = v(q)
    equations(int).v(v, sol.t, q, params)
end


function residual!(b::AbstractVector{ST}, int::GeometricIntegrator{<:ImplicitEuler}) where {ST}
    # temporary variables
    local y::ST

    # get cache for internal stages
    local v = cache(int, ST).v
    local v̄ = cache(int, ST).v̄

    # compute b = - (v-v)
    b .= v .- v̄
end


# Compute stages of implicit Euler methods.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:ImplicitEuler}) where {ST}
    @assert axes(x) == axes(b)

    # copy previous solution from solstep to cache
    reset!(cache(int, ST), sol...)

    # compute stages from nonlinear solver solution x
    components!(x, sol, params, int)

    # compute residual vector
    residual!(b, int)
end


function update!(sol, params, x::AbstractVector{DT}, int::GeometricIntegrator{<:ImplicitEuler}) where {DT}
    # copy previous solution from solstep to cache
    reset!(cache(int, DT), sol...)

    # compute vector field at internal stages
    components!(x, sol, params, int)

    # compute final update
    sol.q .+= timestep(int) .* cache(int, DT).v
end


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:ImplicitEuler, <:AbstractProblemODE})
    # call nonlinear solver
    solve!(nlsolution(int), (b,x) -> residual!(b, x, sol, params, int), solver(int))

    # print solver status
    # println(status(solver))

    # check if solution contains NaNs or error bounds are violated
    # println(meets_stopping_criteria(status(solver)))

    # compute final update
    update!(sol, params, nlsolution(int), int)
end

function integrate_step!(int::GeometricIntegrator{<:ImplicitEuler, <:AbstractProblemODE})
    integrate_step!(current(solstep(int)), history(solstep(int)), parameters(solstep(int)), int)
end
