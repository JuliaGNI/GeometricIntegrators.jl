"""
Solution step for a [`PODEProblem`](@ref).

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
SolutionStepPODE(t::TT, q::AT, p::AT, internal::IT=NamedTuple())
```

"""
mutable struct SolutionStepPODE{
            DT <: Number,
            TT <: Real,
            AT <: AbstractArray{DT},
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

    v::VT
    v̄::VT

    f::FT
    f̄::FT

    internal::IT

    function SolutionStepPODE(t::TT, q::AT, p::AT, internal::IT = NamedTuple()) where {DT, TT, AT <: AbstractArray{DT}, IT}
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

function current(solstep::SolutionStepPODE)
    (t = solstep.t, q = solstep.q, p = solstep.p)
end

function previous(solstep::SolutionStepPODE)
    (t = solstep.t̄, q = solstep.q̄, p = solstep.p̄)
end

function Base.copy!(solstep::SolutionStepPODE, sol::NamedTuple)
    solstep.t  = sol.t
    solstep.q .= sol.q
    solstep.p .= sol.p
    solstep.q̃ .= 0
    solstep.v .= 0
    solstep.f .= 0
    return solstep
end

function GeometricBase.reset!(solstep::SolutionStepPODE)
    solstep.t̄  = solstep.t
    solstep.q̄ .= solstep.q
    solstep.p̄ .= solstep.p
    solstep.v̄ .= solstep.v
    solstep.f̄ .= solstep.f
    return solstep
end

function update!(solstep::SolutionStepPODE{DT,TT,AT}, Δt::TT, Δq::AT, Δp::AT) where {DT,TT,AT}
    solstep.t += Δt
    for k in eachindex(Δq,Δp)
        solstep.q[k], solstep.q̃[k] = compensated_summation(Δq[k], solstep.q[k], solstep.q̃[k])
        solstep.p[k], solstep.p̃[k] = compensated_summation(Δp[k], solstep.p[k], solstep.p̃[k])
    end
    return solstep
end
