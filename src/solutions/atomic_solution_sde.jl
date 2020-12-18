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
* `ΔW`: Wiener process driving the stochastic process q
* `ΔZ`: Wiener process driving the stochastic process q
* `K`:  integer parameter defining the truncation of the increments of the Wiener process (for strong solutions)
"""
mutable struct AtomicSolutionSDE{DT <: Number, TT <: Real, AT <: AbstractArray{DT}, IT <: NamedTuple} <: AtomicSolution{DT,TT,AT}
    t::TT
    t̅::TT

    q::AT
    q̅::AT
    q̃::AT

    ΔW::AT
    ΔZ::AT

    K::Int

    internal::IT

    function AtomicSolutionSDE{DT,TT,AT,IT}(nd, nm, internal::IT) where {DT,TT,AT,IT}
        new(zero(TT), zero(TT),
            AT(zeros(DT, nd)), AT(zeros(DT, nd)), AT(zeros(DT, nd)),
            AT(zeros(DT, nm)), AT(zeros(DT, nm)), 0,
            internal)
    end

    function AtomicSolutionSDE{DT,TT,AT,IT}(t::TT, q::AT, ΔW::AT, ΔZ::AT, internal::IT) where {DT,TT,AT,IT}
        new(zero(t), zero(t),
            zero(q), zero(q), zero(q),
            zero(ΔW), zero(ΔZ), 0,
            internal)
    end
end

AtomicSolutionSDE(::Type{DT}, ::Type{TT}, ::Type{AT}, nd, nm, internal::IT=NamedTuple()) where {DT,TT,AT,IT} = AtomicSolutionSDE{DT,TT,AT,IT}(nd, nm, internal)
AtomicSolutionSDE(t::TT, q::AT, ΔW::AT, ΔZ::AT, internal::IT=NamedTuple()) where {DT, TT, AT <: AbstractArray{DT}, IT} = AtomicSolutionSDE{DT,TT,AT,IT}(t, q, ΔW, ΔZ, internal)

function set_solution!(asol::AtomicSolutionSDE, sol)
    t, q = sol
    asol.t  = t
    asol.q .= q
end

function set_increments!(asol::AtomicSolutionSDE, incs)
    ΔW, ΔZ = incs
    asol.ΔW .= ΔW
    asol.ΔZ .= ΔZ
end

function get_solution(asol::AtomicSolutionSDE)
    (asol.t, asol.q)
end

function get_increments(asol::AtomicSolutionSDE)
    (asol.ΔW, asol.ΔZ)
end

function get_increments!(asol::AtomicSolutionSDE, ΔW, ΔZ)
    ΔW .= asol.ΔW
    ΔZ .= asol.ΔZ
end

function Common.reset!(asol::AtomicSolutionSDE, Δt)
    asol.t̅  = asol.t
    asol.q̅ .= asol.q
    asol.t += Δt
end

function update!(asol::AtomicSolutionSDE{DT}, y::Vector{DT}) where {DT}
    for k in eachindex(y)
        update!(asol, y[k], k)
    end
end

function update!(asol::AtomicSolutionSDE{DT}, y::DT, k::Union{Int,CartesianIndex}) where {DT}
    asol.q[k], asol.q̃[k] = compensated_summation(y, asol.q[k], asol.q̃[k])
end
