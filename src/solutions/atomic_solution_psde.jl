"""
Atomic solution for an SDE.

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
* `ΔW`: Wiener process driving the stochastic process q
* `ΔZ`: Wiener process driving the stochastic process q
* `K`:  integer parameter defining the truncation of the increments of the Wiener process (for strong solutions)
"""
mutable struct AtomicSolutionPSDE{DT <: Number, TT <: Real, AT <: AbstractArray{DT}, IT <: NamedTuple} <: AtomicSolution{DT,TT,AT}
    t::TT
    t̅::TT

    q::AT
    q̅::AT
    q̃::AT

    p::AT
    p̅::AT
    p̃::AT

    ΔW::AT
    ΔZ::AT

    K::Int

    internal::IT

    function AtomicSolutionPSDE{DT,TT,AT,IT}(nd, nm, internal::IT) where {DT,TT,AT,IT}
        new(zero(TT), zero(TT),
            AT(zeros(DT, nd)), AT(zeros(DT, nd)), AT(zeros(DT, nd)),
            AT(zeros(DT, nd)), AT(zeros(DT, nd)), AT(zeros(DT, nd)),
            AT(zeros(DT, nm)), AT(zeros(DT, nm)), 0, internal)
    end

    function AtomicSolutionPSDE{DT,TT,AT,IT}(t::TT, q::AT, p::AT, ΔW::AT, ΔZ::AT, internal::IT) where {DT,TT,AT,IT}
        new(zero(t), zero(t),
            zero(q), zero(q), zero(q),
            zero(p), zero(p), zero(p),
            zero(ΔW), zero(ΔZ), 0,
            internal)
    end
end

AtomicSolutionPSDE(::Type{DT}, ::Type{TT}, ::Type{AT}, nd, nm, internal::IT=NamedTuple()) where {DT,TT,AT,IT} = AtomicSolutionPSDE{DT,TT,AT,IT}(nd, nm, internal)
AtomicSolutionPSDE(t::TT, q::AT, p::AT, ΔW::AT, ΔZ::AT, internal::IT=NamedTuple()) where {DT, TT, AT <: AbstractArray{DT}, IT} = AtomicSolutionPSDE{DT,TT,AT,IT}(t, q, p, ΔW, ΔZ, internal)

function set_solution!(asol::AtomicSolutionPSDE, sol)
    t, q, p = sol
    asol.t  = t
    asol.q .= q
    asol.p .= p
end

function set_increments!(asol::AtomicSolutionPSDE, incs)
    ΔW, ΔZ = incs
    asol.ΔW .= ΔW
    asol.ΔZ .= ΔZ
end

function get_solution(asol::AtomicSolutionPSDE)
    (asol.t, asol.q, asol.p)
end

function get_increments(asol::AtomicSolutionPSDE)
    (asol.ΔW, asol.ΔZ)
end

function get_increments!(asol::AtomicSolutionPSDE, ΔW, ΔZ)
    ΔW .= asol.ΔW
    ΔZ .= asol.ΔZ
end

function Common.reset!(asol::AtomicSolutionPSDE, Δt)
    asol.t̅  = asol.t
    asol.q̅ .= asol.q
    asol.p̅ .= asol.p
    asol.t += Δt
end

function update!(asol::AtomicSolutionPSDE{DT}, y::Vector{DT}, z::Vector{DT}) where {DT}
    for k in eachindex(y,z)
        update!(asol, y[k], z[k], k)
    end
end

function update!(asol::AtomicSolutionPSDE{DT}, y::DT, z::DT, k::Union{Int,CartesianIndex}) where {DT}
    asol.q[k], asol.q̃[k] = compensated_summation(y, asol.q[k], asol.q̃[k])
    asol.p[k], asol.p̃[k] = compensated_summation(z, asol.p[k], asol.p̃[k])
end
