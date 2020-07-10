"""
Atomic solution for an ODE.

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
* `v`: vector field of q
* `v̅`: vector field of q̅
"""
mutable struct AtomicSolutionODE{DT <: Number, TT <: Real, AT <: AbstractArray{DT}, IT <: NamedTuple} <: AtomicSolution{DT,TT}
    t::TT
    t̅::TT

    q::AT
    q̅::AT
    q̃::AT

    v::AT
    v̅::AT

    internal::IT

    function AtomicSolutionODE{DT,TT,AT,IT}(nd, internal::IT) where {DT,TT,AT,IT}
        new(zero(TT), zero(TT),
            AT(zeros(DT, nd)), AT(zeros(DT, nd)), AT(zeros(DT, nd)),
            AT(zeros(DT, nd)), AT(zeros(DT, nd)),
            internal)
    end

    function AtomicSolutionODE{DT,TT,AT,IT}(t::TT, q::AT, internal::IT) where {DT,TT,AT,IT}
        new(zero(t), zero(t),
            zero(q), zero(q), zero(q),
            zero(q), zero(q),
            internal)
    end
end

AtomicSolutionODE(DT, TT, AT, nd, internal::IT=NamedTuple()) where {IT} = AtomicSolutionODE{DT,TT,AT,IT}(nd, internal)
AtomicSolutionODE(t::TT, q::AT, internal::IT=NamedTuple()) where {DT, TT, AT <: AbstractArray{DT}, IT} = AtomicSolutionODE{DT,TT,AT,IT}(t, q, internal)


function set_initial_conditions!(asol::AtomicSolutionODE, equ::AbstractEquationODE)
    t, q = initial_conditions(equ)
    asol.t  = t
    asol.q .= q
    asol.v .= 0
end

function set_solution!(asol::AtomicSolutionODE, sol)
    t, q = sol
    asol.t  = t
    asol.q .= q
    asol.v .= 0
end

function get_solution(asol::AtomicSolutionODE)
    (asol.t, asol.q)
end

function Common.reset!(asol::AtomicSolutionODE, Δt)
    asol.t̅  = asol.t
    asol.q̅ .= asol.q
    asol.v̅ .= asol.v
    asol.t += Δt
end

function update!(asol::AtomicSolutionODE{DT}, y::Vector{DT}) where {DT}
    for k in eachindex(y)
        update!(asol, y[k], k)
    end
end

function update!(asol::AtomicSolutionODE{DT}, y::DT, k::Union{Int,CartesianIndex}) where {DT}
    asol.q[k], asol.q̃[k] = compensated_summation(y, asol.q[k], asol.q̃[k])
end
