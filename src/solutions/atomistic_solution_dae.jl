"""
Atomistic solution for an DAE.

### Fields

* `t`: time of current time step
* `t̅`: time of previous time step
* `q`: current solution of q
* `q̅`: previous solution of q
* `q̃`: compensated summation error of q
* `λ`: current solution of λ
* `λ̅`: previous solution of λ
* `v`: vector field of q
* `v̅`: vector field of q̅
* `u`: projective vector field of q
* `u̅`: projective vector field of q̅
"""
mutable struct AtomisticSolutionDAE{DT,TT} <: AtomisticSolution{DT,TT}
    t::TT
    t̅::TT

    q::Vector{DT}
    q̅::Vector{DT}
    q̃::Vector{DT}

    λ::Vector{DT}
    λ̅::Vector{DT}

    v::Vector{DT}
    v̅::Vector{DT}

    u::Vector{DT}
    u̅::Vector{DT}

    function AtomisticSolutionDAE{DT, TT}(nd) where {DT <: Number, TT <: Real}
        new(zero(TT), zero(TT), zeros(DT, nd), zeros(DT, nd), zeros(DT, nd),
                                zeros(DT, nd), zeros(DT, nd),
                                zeros(DT, nd), zeros(DT, nd),
                                zeros(DT, nd), zeros(DT, nd))
    end
end

AtomisticSolutionDAE(DT, TT, nd) = AtomisticSolutionDAE{DT, TT}(nd)

function CommonFunctions.reset!(cache::AtomisticSolutionDAE, Δt)
    cache.t̅  = cache.t
    cache.q̅ .= cache.q
    cache.λ̅ .= cache.λ
    cache.v̅ .= cache.v
    cache.u̅ .= cache.u
    cache.t += Δt
end

function update!(asol::AtomisticSolutionDAE{DT}, v::Vector{DT}, λ::Vector{DT}) where {DT}
    for k in eachindex(v)
        update!(asol, v[k], λ[k])
    end
end

function update!(asol::AtomisticSolutionDAE{DT}, v::DT, λ::DT, k::Int) where {DT}
    asol.q[k], asol.q̃[k] = compensated_summation(v, asol.q[k], asol.q̃[k])
    asol.λ[k] = λ
end
