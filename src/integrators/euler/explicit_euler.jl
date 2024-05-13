"""
Explicit Euler Method.

$(reference(Val(:ExplicitEuler)))
"""
struct ExplicitEuler <: ODEMethod end

isexplicit(method::ExplicitEuler) = true
isimplicit(method::ExplicitEuler) = false
issymmetric(method::ExplicitEuler) = false
issymplectic(method::ExplicitEuler) = false

@doc raw"""
Implicit Euler integrator cache.
"""
struct ExplicitEulerCache{DT,D} <: ODEIntegratorCache{DT,D}
    v::Vector{DT}

    function ExplicitEulerCache{DT,D}() where {DT,D}
        v = zeros(DT, D)
        new(v)
    end
end

function Cache{ST}(problem::AbstractProblem, method::ExplicitEuler; kwargs...) where {ST}
    ExplicitEulerCache{ST, ndims(problem)}(; kwargs...)
end

@inline CacheType(ST, problem::AbstractProblem, method::ExplicitEuler) = ExplicitEulerCache{ST, ndims(problem)}


function update!(sol, params, _, int::GeometricIntegrator{<:ExplicitEuler})
    # compute final update
    sol.q .+= timestep(int) .* cache(int).v
end


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:ExplicitEuler, <:AbstractProblemODE})
    # compute vector field
    equations(int).v(cache(int).v, sol.t, sol.q, params)

    # compute final update
    update!(sol, params, nothing, int)
end

function integrate_step!(int::GeometricIntegrator{<:ExplicitEuler, <:AbstractProblemODE})
    integrate_step!(current(solstep(int)), history(solstep(int)), parameters(solstep(int)), int)
end
