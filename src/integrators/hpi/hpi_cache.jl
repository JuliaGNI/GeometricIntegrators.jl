@doc raw"""
Cache for variational integrator in position-momentum form.

### Fields

* `q`: internal stages of solution
* `Θ`: implicit function evaluated on solution
* `f`: vector field of implicit function
"""
struct HPICache{DT,D,A} <: IODEIntegratorCache{DT,D}
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

    function HPICache{DT,D,A}() where {DT,D,A}
        new(zeros(DT, D+A), zeros(DT, A),
            zeros(DT,D), zeros(DT,D), zeros(DT,D), zeros(DT,D),
            zeros(DT,D), zeros(DT,D), zeros(DT,D), zeros(DT,D),
            zeros(DT,D), zeros(DT,D), zeros(DT,D), zeros(DT,D),
            zeros(DT,D,D), zeros(DT,D,D), zeros(DT,A,D)
        )
    end
end

nlsolution(cache::HPICache) = cache.x

function Cache{ST}(problem::AbstractProblemIODE, method::HPIMethod; kwargs...) where {ST}
    HPICache{ST, ndims(problem), nparams(method)}(; kwargs...)
end

@inline CacheType(ST, problem::AbstractProblemIODE, method::HPIMethod) = HPICache{ST, ndims(problem), nparams(method)}
