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
AtomicSolutionODE(t::TT, q::AT, internal::IT=NamedTuple())
```

"""
mutable struct AtomicSolutionODE{
            DT <: Number, 
            TT <: Real, 
            AT <: AbstractArray{DT}, 
            VT <: AbstractArray{DT}, 
            IT <: NamedTuple} <: AtomicSolution{DT,TT,AT}

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
        new{DT,TT,AT,typeof(v),IT}(
            copy(t), zero(t),
            copy(q), zero(q), zero(q), 
            v, v̄, internal)
    end
end

function current(asol::AtomicSolutionODE)
    (t = asol.t, q = asol.q)
end

function previous(asol::AtomicSolutionODE)
    (t = asol.t̄, q = asol.q̄)
end

function Base.copy!(asol::AtomicSolutionODE, sol::NamedTuple)
    asol.t  = sol.t
    asol.t̄  = 0
    asol.q .= sol.q
    asol.q̄ .= 0
    asol.q̃ .= 0
    asol.v .= 0
    asol.v̄ .= 0
    return asol
end

function GeometricBase.reset!(asol::AtomicSolutionODE)
    asol.t̄  = asol.t
    asol.q̄ .= asol.q
    asol.v̄ .= asol.v
    return asol
end

function update!(asol::AtomicSolutionODE{DT,TT,AT}, Δt::TT, Δq::AT) where {DT,TT,AT}
    asol.t += Δt
    for k in eachindex(Δq)
        asol.q[k], asol.q̃[k] = compensated_summation(Δq[k], asol.q[k], asol.q̃[k])
    end
    return asol
end
