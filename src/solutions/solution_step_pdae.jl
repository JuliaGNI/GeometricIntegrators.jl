"""
Solution step for a [`PDAEProblem`](@ref).

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
SolutionStepPDAE(t::TT, q::AT, p::AT, λ::AT, internal::IT=NamedTuple())
```

"""
mutable struct SolutionStepPDAE{
                DT <: Number,
                TT <: Real,
                AT <: AbstractArray{DT},
                ΛT <: AbstractArray{DT},
                VT <: AbstractArray{DT},
                FT <: AbstractArray{DT},
                IT <: NamedTuple} <: SolutionStep{DT,TT,AT}
    
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

    function SolutionStepPDAE(t::TT, q::AT, p::AT, λ::ΛT, internal::IT = NamedTuple()) where {DT, TT, AT <: AbstractArray{DT}, ΛT <: AbstractArray{DT}, IT}
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

function current(solstep::SolutionStepPDAE)
    (t = solstep.t, q = solstep.q, p = solstep.p, λ = solstep.λ)
end

function previous(solstep::SolutionStepPDAE)
    (t = solstep.t̄, q = solstep.q̄, p = solstep.p̄, λ = solstep.λ̄)
end

function Base.copy!(solstep::SolutionStepPDAE, sol::NamedTuple)
    solstep.t  = sol.t
    solstep.q .= sol.q
    solstep.p .= sol.p
    solstep.λ .= sol.λ
    solstep.q̃ .= 0
    solstep.v .= 0
    solstep.f .= 0
    solstep.u .= 0
    solstep.g .= 0
    return solstep
end

function GeometricBase.reset!(solstep::SolutionStepPDAE)
    solstep.t̄  = solstep.t
    solstep.q̄ .= solstep.q
    solstep.p̄ .= solstep.p
    solstep.λ̄ .= solstep.λ
    solstep.v̄ .= solstep.v
    solstep.f̄ .= solstep.f
    solstep.ū .= solstep.u
    solstep.ḡ .= solstep.g
    return solstep
end

function update!(solstep::SolutionStepPDAE{DT,TT,AT,ΛT}, Δt::TT, Δq::AT, Δp::AT, λ::ΛT) where {DT,TT,AT,ΛT}
    solstep.t += Δt
    for k in eachindex(Δq,Δp)
        solstep.q[k], solstep.q̃[k] = compensated_summation(Δq[k], solstep.q[k], solstep.q̃[k])
        solstep.p[k], solstep.p̃[k] = compensated_summation(Δp[k], solstep.p[k], solstep.p̃[k])
    end
    solstep.λ .= λ
    return solstep
end
