"""
Atomic solution for an PODE.

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
* `p`: current solution of p
* `p̄`: previous solution of p
* `p̃`: compensated summation error of p
* `v`: vector field of q
* `v̄`: vector field of q̄
* `f`: vector field of p
* `f̄`: vector field of p̄
* `internal`: internal variables of the integrator (e.g., internal stages of a Runge-Kutta methods or solver output)

### Constructors

```julia
AtomicSolutionPODE(t::TT, q::AT, p::AT, internal::IT=NamedTuple())
```

"""
mutable struct AtomicSolutionPODE{
            DT <: Number,
            TT <: Real,
            AT <: AbstractArray{DT},
            VT <: AbstractArray{DT},
            FT <: AbstractArray{DT},
            IT <: NamedTuple} <: AtomicSolution{DT,TT,AT}

    t::TT
    t̄::TT

    q::AT
    q̄::AT
    q̃::AT

    p::AT
    p̄::AT
    p̃::AT

    v::VT
    v̄::VT

    f::FT
    f̄::FT

    internal::IT

    function AtomicSolutionPODE(t::TT, q::AT, p::AT, internal::IT = NamedTuple()) where {DT, TT, AT <: AbstractArray{DT}, IT}
        v = vectorfield(q)
        v̄ = vectorfield(q)
        f = vectorfield(p)
        f̄ = vectorfield(p)

        new{DT,TT,AT,typeof(v),typeof(f),IT}(
            copy(t), zero(t),
            copy(q), zero(q), zero(q),
            copy(p), zero(p), zero(p),
            v, v̄, f, f̄, internal
        )
    end
end

function current(asol::AtomicSolutionPODE)
    (t = asol.t, q = asol.q, p = asol.p)
end

function previous(asol::AtomicSolutionPODE)
    (t = asol.t̄, q = asol.q̄, p = asol.p̄)
end

function Base.copy!(asol::AtomicSolutionPODE, sol::NamedTuple)
    asol.t  = sol.t
    asol.q .= sol.q
    asol.p .= sol.p
    asol.q̃ .= 0
    asol.v .= 0
    asol.f .= 0
    return asol
end

function GeometricBase.reset!(asol::AtomicSolutionPODE)
    asol.t̄  = asol.t
    asol.q̄ .= asol.q
    asol.p̄ .= asol.p
    asol.v̄ .= asol.v
    asol.f̄ .= asol.f
    return asol
end

function update!(asol::AtomicSolutionPODE{DT,TT,AT}, Δt::TT, Δq::AT, Δp::AT) where {DT,TT,AT}
    asol.t += Δt
    for k in eachindex(Δq,Δp)
        asol.q[k], asol.q̃[k] = compensated_summation(Δq[k], asol.q[k], asol.q̃[k])
        asol.p[k], asol.p̃[k] = compensated_summation(Δp[k], asol.p[k], asol.p̃[k])
    end
    return asol
end
