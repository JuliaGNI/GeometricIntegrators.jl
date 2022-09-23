"""
Atomic solution for an PDAE.

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
* `λ`: current solution of λ
* `λ̄`: previous solution of λ
* `v`: vector field of q
* `v̄`: vector field of q̄
* `f`: vector field of p
* `f̄`: vector field of p̄
* `u`: projective vector field of q
* `ū`: projective vector field of q̄
* `g`: projective vector field of p
* `ḡ`: projective vector field of p̄
* `internal`: internal variables of the integrator (e.g., internal stages of a Runge-Kutta methods or solver output)

### Constructors

```julia
AtomicSolutionPDAE(t::TT, q::AT, p::AT, λ::AT, internal::IT=NamedTuple())
```

"""
mutable struct AtomicSolutionPDAE{
                DT <: Number,
                TT <: Real,
                AT <: AbstractArray{DT},
                ΛT <: AbstractArray{DT},
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

    λ::ΛT
    λ̄::ΛT

    v::VT
    v̄::VT
    f::FT
    f̄::FT
    
    u::VT
    ū::VT
    g::FT
    ḡ::FT

    internal::IT

    function AtomicSolutionPDAE(t::TT, q::AT, p::AT, λ::ΛT, internal::IT = NamedTuple()) where {DT, TT, AT <: AbstractArray{DT}, ΛT <: AbstractArray{DT}, IT}
        v = vectorfield(q)
        v̄ = vectorfield(q)
        f = vectorfield(p)
        f̄ = vectorfield(p)
        u = vectorfield(q)
        ū = vectorfield(q)
        g = vectorfield(p)
        ḡ = vectorfield(p)

        new{DT,TT,AT,ΛT,typeof(v),typeof(f),IT}(
            copy(t), zero(t),
            copy(q), zero(q), zero(q),
            copy(p), zero(p), zero(p),
            copy(λ), zero(λ),
            v, v̄, f, f̄, u, ū, g, ḡ, internal)
    end
end

function current(asol::AtomicSolutionPDAE)
    (t = asol.t, q = asol.q, p = asol.p, λ = asol.λ)
end

function previous(asol::AtomicSolutionPDAE)
    (t = asol.t̄, q = asol.q̄, p = asol.p̄, λ = asol.λ̄)
end

function Base.copy!(asol::AtomicSolutionPDAE, sol::NamedTuple)
    asol.t  = sol.t
    asol.q .= sol.q
    asol.p .= sol.p
    asol.λ .= sol.λ
    asol.q̃ .= 0
    asol.v .= 0
    asol.f .= 0
    asol.u .= 0
    asol.g .= 0
    return asol
end

function GeometricBase.reset!(asol::AtomicSolutionPDAE)
    asol.t̄  = asol.t
    asol.q̄ .= asol.q
    asol.p̄ .= asol.p
    asol.λ̄ .= asol.λ
    asol.v̄ .= asol.v
    asol.f̄ .= asol.f
    asol.ū .= asol.u
    asol.ḡ .= asol.g
    return asol
end

function update!(asol::AtomicSolutionPDAE{DT,TT,AT,ΛT}, Δt::TT, Δq::AT, Δp::AT, λ::ΛT) where {DT,TT,AT,ΛT}
    asol.t += Δt
    for k in eachindex(Δq,Δp)
        asol.q[k], asol.q̃[k] = compensated_summation(Δq[k], asol.q[k], asol.q̃[k])
        asol.p[k], asol.p̃[k] = compensated_summation(Δp[k], asol.p[k], asol.p̃[k])
    end
    asol.λ .= λ
    return asol
end
