@doc raw"""
Cache for variational integrator in position-momentum form.

### Fields

* `q`: internal stages of solution
* `Θ`: implicit function evaluated on solution
* `f`: vector field of implicit function
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

function reset!(cache::IntegratorCachePMVI, t, q, p)
    copyto!(cache.q̄, q)
    copyto!(cache.p̄, p)
end

nlsolution(cache::IntegratorCachePMVI) = cache.x
