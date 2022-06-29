"""
Atomic solution for an ODE.

### Parameters

* `DT`: data type
* `TT`: time step type
* `AT`: array type
* `IT`: internal variable types

### Fields

* `t`: time of current time step
* `t̄`: time of previous time step
* `q`: current solution of q
* `q̄`: previous solution of q
* `q̃`: compensated summation error of q
* `v`: vector field of q
* `v̄`: vector field of q̄
* `internal`: internal variables of the integrator (e.g., internal stages of a Runge-Kutta methods or solver output)

### Constructors

```julia
AtomicSolutionODE{DT,TT,AT,IT}(nd, internal::IT)
AtomicSolutionODE{DT,TT,AT,IT}(t::TT, q::AT, internal::IT)
AtomicSolutionODE(DT, TT, AT, nd, internal::IT=NamedTuple())
AtomicSolutionODE(t::TT, q::AT, internal::IT=NamedTuple())
```

* `nd`: dimension of the state vector

"""
mutable struct AtomicSolutionODE{DT <: Number, TT <: Real, AT <: AbstractArray{DT}, VT <: AbstractArray{DT}, IT <: NamedTuple} <: AtomicSolution{DT,TT,AT}
    t::TT
    t̄::TT

    q::AT
    q̄::AT
    q̃::AT

    v::VT
    v̄::VT

    internal::IT

    function AtomicSolutionODE(t::TT, q::AT, internal::IT = NamedTuple()) where {DT, TT, AT <: AbstractArray{DT}, IT}
        v = vectorfield(q)
        v̄ = vectorfield(q)
        new{DT,TT,AT,typeof(v),IT}(zero(t), zero(t), zero(q), zero(q), zero(q), v, v̄, internal)
    end
end

AtomicSolutionODE{DT,TT,AT}(nd, args...) where {DT,TT,AT} = AtomicSolutionODE(zero(TT), AT(zeros(DT, nd)), args...)

AtomicSolutionODE(DT, TT, AT, nd, args...) = AtomicSolutionODE{DT,TT,AT}(nd, args...)


function set_initial_conditions!(asol::AtomicSolutionODE, equ::AbstractProblemODE, i::Int=1)
    @assert i ≥ nsamples(equ)
    t, q = initial_conditions(equ)
    asol.t  = t
    asol.q .= q[i]
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

function GeometricBase.reset!(asol::AtomicSolutionODE, Δt)
    asol.t̄  = asol.t
    asol.q̄ .= asol.q
    asol.v̄ .= asol.v
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
