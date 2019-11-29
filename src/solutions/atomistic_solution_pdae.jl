"""
Atomistic solution for an PDAE.

### Fields

* `t`: time of current time step
* `t̅`: time of previous time step
* `q`: current solution of q
* `q̅`: previous solution of q
* `q̃`: compensated summation error of q
* `p`: current solution of p
* `p̅`: previous solution of p
* `p̃`: compensated summation error of p
* `λ`: current solution of λ
* `λ̅`: previous solution of λ
* `v`: vector field of q
* `v̅`: vector field of q̅
* `f`: vector field of p
* `f̅`: vector field of p̅
* `u`: projective vector field of q
* `u̅`: projective vector field of q̅
* `g`: projective vector field of p
* `g̅`: projective vector field of p̅
"""
mutable struct AtomisticSolutionPDAE{DT,TT} <: AtomisticSolution{DT,TT}
    t::TT
    t̅::TT

    q::Vector{DT}
    q̅::Vector{DT}
    q̃::Vector{DT}

    p::Vector{DT}
    p̅::Vector{DT}
    p̃::Vector{DT}

    λ::Vector{DT}
    λ̅::Vector{DT}

    v::Vector{DT}
    v̅::Vector{DT}
    f::Vector{DT}
    f̅::Vector{DT}

    u::Vector{DT}
    u̅::Vector{DT}
    g::Vector{DT}
    g̅::Vector{DT}

    function AtomisticSolutionPDAE{DT, TT}(nd) where {DT <: Number, TT <: Real}
        new(zero(TT), zero(TT), zeros(DT, nd), zeros(DT, nd), zeros(DT, nd),
                                zeros(DT, nd), zeros(DT, nd), zeros(DT, nd),
                                zeros(DT, nd), zeros(DT, nd),
                                zeros(DT, nd), zeros(DT, nd), zeros(DT, nd), zeros(DT, nd),
                                zeros(DT, nd), zeros(DT, nd), zeros(DT, nd), zeros(DT, nd))
    end
end

AtomisticSolutionPDAE(DT, TT, nd) = AtomisticSolutionPDAE{DT, TT}(nd)

function CommonFunctions.reset!(cache::AtomisticSolutionPDAE, Δt)
    cache.t̅  = cache.t
    cache.q̅ .= cache.q
    cache.p̅ .= cache.p
    cache.λ̅ .= cache.λ
    cache.v̅ .= cache.v
    cache.f̅ .= cache.f
    cache.u̅ .= cache.u
    cache.g̅ .= cache.g
    cache.t += Δt
end

function update!(asol::AtomisticSolutionPODE{DT}, v::Vector{DT}, f::Vector{DT}, λ::Vector{DT}) where {DT}
    for k in eachindex(v,f)
        update!(asol, v[k], f[k], λ[k])
    end
end

function update!(asol::AtomisticSolutionPODE{DT}, v::DT, f::DT, λ::DT, k::Int) where {DT}
    asol.q[k], asol.q̃[k] = compensated_summation(v, asol.q[k], asol.q̃[k])
    asol.p[k], asol.p̃[k] = compensated_summation(f, asol.p[k], asol.p̃[k])
    asol.λ[k] = λ
end
