"""
Atomic solution for an DAE.

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
* `λ`: current solution of λ
* `λ̄`: previous solution of λ
* `v`: vector field of q
* `v̄`: vector field of q̄
* `u`: projective vector field of q
* `ū`: projective vector field of q̄
* `internal`: internal variables of the integrator (e.g., internal stages of a Runge-Kutta methods or solver output)

### Constructors

```julia
AtomicSolutionDAE(t::TT, q::AT, λ::AT, internal::IT=NamedTuple())
```

"""
mutable struct AtomicSolutionDAE{
            DT <: Number,
            TT <: Real,
            AT <: AbstractArray{DT},
            ΛT <: AbstractArray{DT},
            VT <: AbstractArray{DT},
            IT <: NamedTuple} <: AtomicSolution{DT,TT,AT}

    t::TT
    t̄::TT

    q::AT
    q̄::AT
    q̃::AT

    λ::ΛT
    λ̄::ΛT

    v::VT
    v̄::VT

    u::VT
    ū::VT

    internal::IT

    function AtomicSolutionDAE(t::TT, q::AT, λ::ΛT, internal::IT = NamedTuple()) where {DT, TT, AT <: AbstractArray{DT}, ΛT <: AbstractArray{DT}, IT}
        v = vectorfield(q)
        v̄ = vectorfield(q)
        u = vectorfield(q)
        ū = vectorfield(q)

        new{DT,TT,AT,ΛT,typeof(v),IT}(
            copy(t), zero(t),
            copy(q), zero(q), zero(q),
            copy(λ), zero(λ),
            v, v̄, u, ū, internal)
    end
end

function current(asol::AtomicSolutionDAE)
    (t = asol.t, q = asol.q, λ = asol.λ)
end

function previous(asol::AtomicSolutionDAE)
    (t = asol.t̄, q = asol.q̄, λ = asol.λ̄)
end

function Base.copy!(asol::AtomicSolutionDAE, sol::NamedTuple)
    asol.t  = sol.t
    asol.q .= sol.q
    asol.λ .= sol.λ
    asol.q̃ .= 0
    asol.v .= 0
    asol.u .= 0
    return asol
end

function GeometricBase.reset!(asol::AtomicSolutionDAE)
    asol.t̄  = asol.t
    asol.q̄ .= asol.q
    asol.λ̄ .= asol.λ
    asol.v̄ .= asol.v
    asol.ū .= asol.u
    return asol
end

function update!(asol::AtomicSolutionDAE{DT}, Δt, Δq::Vector{DT}, λ::Vector{DT}) where {DT}
    asol.t += Δt
    for k in eachindex(Δq)
        asol.q[k], asol.q̃[k] = compensated_summation(Δq[k], asol.q[k], asol.q̃[k])
    end
    asol.λ .= λ
    return asol
end
