@doc raw"""
Cache for variational integrator in position-momentum form.

### Fields

* `x`: nonlinear solver solution vector
* `q`: solution at next timestep
* `p`: momentum at next timestep
* `θ`: implicit function evaluated on solution at next timestep
* `f`: vector field of implicit function evaluated on solution at next timestep
* `q̄`: solution at previous timestep
* `p̄`: momentum at previous timestep
* `θ̄`: implicit function evaluated on solution at previous timestep
* `f̄`: vector field of implicit function evaluated on solution at previous timestep
"""
struct IntegratorCachePMVI{DT,D} <: IODEIntegratorCache{DT,D}
    x::Vector{DT}

    q::Vector{DT}
    p::Vector{DT}
    θ::Vector{DT}
    f::Vector{DT}

    q̄::Vector{DT}
    p̄::Vector{DT}
    θ̄::Vector{DT}
    f̄::Vector{DT}

    q̃::Vector{DT}
    ṽ::Vector{DT}
    θ̃::Vector{DT}
    f̃::Vector{DT}

    function IntegratorCachePMVI{DT,D}() where {DT,D}
        new(zeros(DT,D), 
            zeros(DT,D), zeros(DT,D), zeros(DT,D), zeros(DT,D),
            zeros(DT,D), zeros(DT,D), zeros(DT,D), zeros(DT,D),
            zeros(DT,D), zeros(DT,D), zeros(DT,D), zeros(DT,D))
    end
end

nlsolution(cache::IntegratorCachePMVI) = cache.x

function Cache{ST}(problem::AbstractProblemIODE, method::PMVIMethod; kwargs...) where {ST}
    IntegratorCachePMVI{ST, ndims(problem)}(; kwargs...)
end

@inline CacheType(ST, problem::AbstractProblemIODE, ::PMVIMethod) = IntegratorCachePMVI{ST, ndims(problem)}
