"""
Atomic solution for an PDAE.

### Parameters

* `DT`: data type
* `TT`: time step type
* `AT`: array type
* `IT`: internal variable types

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
mutable struct AtomicSolutionPDAE{DT, TT, AT <: AbstractArray{DT}, IT <: NamedTuple} <: AtomicSolution{DT,TT}
    t::TT
    t̅::TT

    q::AT
    q̅::AT
    q̃::AT

    p::AT
    p̅::AT
    p̃::AT

    λ::AT
    λ̅::AT

    v::AT
    v̅::AT
    f::AT
    f̅::AT

    u::AT
    u̅::AT
    g::AT
    g̅::AT

    internal::IT

    function AtomicSolutionPDAE{DT,TT,AT,IT}(nd, nm, internal::IT) where {DT <: Number, TT <: Real, AT, IT <: NamedTuple}
        new(zero(TT), zero(TT),
            zeros(DT, nd), zeros(DT, nd), zeros(DT, nd),
            zeros(DT, nd), zeros(DT, nd), zeros(DT, nd),
            zeros(DT, nm), zeros(DT, nm),
            zeros(DT, nd), zeros(DT, nd), zeros(DT, nd), zeros(DT, nd),
            zeros(DT, nd), zeros(DT, nd), zeros(DT, nd), zeros(DT, nd),
            internal)
    end

    function AtomicSolutionPDAE{DT,TT,AT,IT}(t::TT, q::AT, p::AT, λ::AT, internal::IT) where {DT <: Number, TT <: Real, AT <: AbstractArray{DT}, IT <: NamedTuple}
        new(zero(t), zero(t),
            zero(q), zero(q), zero(q),
            zero(p), zero(p), zero(p),
            zero(λ), zero(λ),
            zero(q), zero(q), zero(p), zero(p),
            zero(λ), zero(λ), zero(λ), zero(λ),
            internal)
    end
end

AtomicSolutionPDAE(DT, TT, AT, nd, nm, internal::IT=NamedTuple()) where {IT} = AtomicSolutionPDAE{DT,TT,AT,IT}(nd, nm, internal)
AtomicSolutionPDAE(t::TT, q::AT, p::AT, λ::AT, internal::IT=NamedTuple()) where {DT, TT, AT <: AbstractArray{DT}, IT} = AtomicSolutionPDAE{DT,TT,AT,IT}(t, q, p, λ, internal)

function set_solution!(asol::AtomicSolutionPDAE, sol)
    t, q, p, λ = sol
    asol.t  = t
    asol.q .= q
    asol.p .= p
    asol.λ .= λ
    asol.v .= 0
    asol.f .= 0
end

function get_solution(asol::AtomicSolutionPDAE)
    (asol.t, asol.q, asol.p, asol.λ)
end

function Common.reset!(asol::AtomicSolutionPDAE, Δt)
    asol.t̅  = asol.t
    asol.q̅ .= asol.q
    asol.p̅ .= asol.p
    asol.λ̅ .= asol.λ
    asol.v̅ .= asol.v
    asol.f̅ .= asol.f
    asol.u̅ .= asol.u
    asol.g̅ .= asol.g
    asol.t += Δt
end

function update!(asol::AtomicSolutionPDAE{DT}, y::Vector{DT}, z::Vector{DT}, λ::Vector{DT}) where {DT}
    for k in eachindex(y,z)
        update!(asol, y[k], z[k], k)
    end
    for k in eachindex(λ)
        asol.λ[k] = λ[k]
    end
end

function update!(asol::AtomicSolutionPDAE{DT}, y::DT, z::DT, k::Union{Int,CartesianIndex}) where {DT}
    asol.q[k], asol.q̃[k] = compensated_summation(y, asol.q[k], asol.q̃[k])
    asol.p[k], asol.p̃[k] = compensated_summation(z, asol.p[k], asol.p̃[k])
end
