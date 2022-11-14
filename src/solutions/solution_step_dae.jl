"""
Solution step for a [`DAEProblem`](@ref).

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
SolutionStepDAE(t::TT, q::AT, λ::AT, internal::IT=NamedTuple())
```

"""
mutable struct SolutionStepDAE{
            DT <: Number,
            TT <: Real,
            AT <: AbstractArray{DT},
            ΛT <: AbstractArray{DT},
            VT <: AbstractArray{DT},
            IT <: NamedTuple} <: SolutionStep{DT,TT,AT}

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

    function SolutionStepDAE(t::TT, q::AT, λ::ΛT, internal::IT = NamedTuple()) where {DT, TT, AT <: AbstractArray{DT}, ΛT <: AbstractArray{DT}, IT}
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

function current(solstep::SolutionStepDAE)
    (t = solstep.t, q = solstep.q, λ = solstep.λ)
end

function previous(solstep::SolutionStepDAE)
    (t = solstep.t̄, q = solstep.q̄, λ = solstep.λ̄)
end

function Base.copy!(solstep::SolutionStepDAE, sol::NamedTuple)
    solstep.t  = sol.t
    solstep.q .= sol.q
    solstep.λ .= sol.λ
    solstep.q̃ .= 0
    solstep.v .= 0
    solstep.u .= 0
    return solstep
end

function GeometricBase.reset!(solstep::SolutionStepDAE)
    solstep.t̄  = solstep.t
    solstep.q̄ .= solstep.q
    solstep.λ̄ .= solstep.λ
    solstep.v̄ .= solstep.v
    solstep.ū .= solstep.u
    return solstep
end

function update!(solstep::SolutionStepDAE{DT}, Δt, Δq::Vector{DT}, λ::Vector{DT}) where {DT}
    solstep.t += Δt
    for k in eachindex(Δq)
        solstep.q[k], solstep.q̃[k] = compensated_summation(Δq[k], solstep.q[k], solstep.q̃[k])
    end
    solstep.λ .= λ
    return solstep
end
