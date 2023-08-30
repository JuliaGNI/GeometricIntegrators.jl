@doc raw"""
Cache for variational integrator in position-momentum form.

### Fields

* `q`: internal stages of solution
* `Θ`: implicit function evaluated on solution
* `f`: vector field of implicit function
"""
struct IntegratorCacheHPI{DT,D,A} <: IODEIntegratorCache{DT,D}
    x::Vector{DT}
    a::Vector{DT}

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

    D₁ϕ::Matrix{DT}
    D₂ϕ::Matrix{DT}
    Dₐϕ::Matrix{DT}

    function IntegratorCacheHPI{DT,D,A}() where {DT,D,A}
        new(zeros(DT, D+A), zeros(DT, A),
            zeros(DT,D), zeros(DT,D), zeros(DT,D), zeros(DT,D),
            zeros(DT,D), zeros(DT,D), zeros(DT,D), zeros(DT,D),
            zeros(DT,D), zeros(DT,D), zeros(DT,D), zeros(DT,D),
            zeros(DT,D,D), zeros(DT,D,D), zeros(DT,A,D)
        )
    end
end

function reset!(cache::IntegratorCacheHPI, t, q, p)
    copyto!(cache.q̄, q)
    copyto!(cache.p̄, p)
end

nlsolution(cache::IntegratorCacheHPI) = cache.x
