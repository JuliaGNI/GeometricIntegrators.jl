@doc raw"""
Degenerate variational integrator cache.

### Fields

* `x`: nonlinear solver solution vector
* `q`: solution at next timestep
* `p`: momentum at next timestep
* `θ`: implicit function evaluated on solution at next timestep
* `v`: vector field evaluated on solution at next timestep
* `f`: force field evaluated on solution at next timestep
* `q̄`: solution at previous timestep
* `p̄`: momentum at previous timestep
* `θ̄`: implicit function evaluated on solution at previous timestep
* `v̄`: vector field at previous timestep
* `f̄`: foce field evaluated on solution at previous timestep
"""
struct DVICache{DT,D} <: IODEIntegratorCache{DT,D}
    x::Vector{DT}

    q::Vector{DT}
    p::Vector{DT}
    v::Vector{DT}
    f::Vector{DT}

    q̄::Vector{DT}
    p̄::Vector{DT}
    v̄::Vector{DT}
    f̄::Vector{DT}

    q̃::Vector{DT}
    p̃::Vector{DT}
    f̃::Vector{DT}

    function DVICache{DT,D}() where {DT,D}
        new(zeros(DT,2D), 
            zeros(DT,D), zeros(DT,D), zeros(DT,D), zeros(DT,D),
            zeros(DT,D), zeros(DT,D), zeros(DT,D), zeros(DT,D),
            zeros(DT,D), zeros(DT,D), zeros(DT,D))
    end
end

nlsolution(cache::DVICache) = cache.x

function Cache{ST}(problem::AbstractProblemIODE, method::DVIMethod; kwargs...) where {ST}
    DVICache{ST, ndims(problem)}(; kwargs...)
end

@inline CacheType(ST, problem::AbstractProblemIODE, method::DVIMethod) = DVICache{ST, ndims(problem)}
