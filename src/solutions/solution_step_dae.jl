"""
Solution step for a [`DAEProblem`](@ref).

### Parameters

* `DT`: data type
* `TT`: time step type
* `AT`: array type
* `IT`: internal variable types

### Fields

* `t`: time of current time step
* `t̄`: time of previous time steps
* `q`: current solution of q
* `q̄`: previous solutions of q
* `q̃`: compensated summation error of q
* `λ`: current solution of λ
* `λ̄`: previous solutions of λ
* `v`: vector field of q
* `v̄`: vector fields of q̄
* `u`: projective vector field of q
* `ū`: projective vector fields of q̄
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
            IT <: NamedTuple,
            NT} <: SolutionStep{DT,TT,AT}

    t::TT
    q::AT
    λ::ΛT
    v::VT
    u::VT

    t̄::OffsetVector{TT, Vector{TT}}
    q̄::OffsetVector{AT, Vector{AT}}
    λ̄::OffsetVector{ΛT, Vector{ΛT}}
    v̄::OffsetVector{VT, Vector{VT}}
    ū::OffsetVector{VT, Vector{VT}}

    q̃::AT

    internal::IT

    function SolutionStepDAE(t::TT, q::AT, λ::ΛT; history = 1, internal::IT = NamedTuple()) where {DT, TT, AT <: AbstractArray{DT}, ΛT <: AbstractArray{DT}, IT}
        @assert history ≥ 1

        t̄ = OffsetVector([zero(t) for _ in 1:history+1], 0:history)
        q̄ = OffsetVector([zero(q) for _ in 1:history+1], 0:history)
        λ̄ = OffsetVector([zero(λ) for _ in 1:history+1], 0:history)
        v̄ = OffsetVector([vectorfield(q) for _ in 1:history+1], 0:history)
        ū = OffsetVector([vectorfield(q) for _ in 1:history+1], 0:history)

        t̄[0]  = t
        q̄[0] .= q
        λ̄[0] .= λ

        new{DT,TT,AT,ΛT,typeof(v̄[0]),IT,history}(t̄[0], q̄[0], λ̄[0], v̄[0], ū[0], t̄, q̄, λ̄ , v̄, ū, zero(q), internal)
    end
end


nhistory(::SolutionStepDAE{DT,TT,AT,ΛT,VT,IT,NT}) where {DT,TT,AT,ΛT,VT,IT,NT} = NT

current(solstep::SolutionStepDAE) = (t = solstep.t, q = solstep.q, λ = solstep.λ)
previous(solstep::SolutionStepDAE) = (t = solstep.t̄[1], q = solstep.q̄[1], λ = solstep.λ̄[1])
history(solstep::SolutionStepDAE) = (t = solstep.t̄, q = solstep.q̄, λ = solstep.λ̄)


function Base.copy!(solstep::SolutionStepDAE, sol::NamedTuple)
    solstep.t  = sol.t
    solstep.q .= sol.q
    solstep.λ .= sol.λ
    solstep.q̃ .= 0
    solstep.v .= 0

    for i in eachhistory(solstep)
        solstep.t̄[i]  = 0
        solstep.q̄[i] .= 0
        solstep.λ̄[i] .= 0
        solstep.v̄[i] .= 0
        solstep.ū[i] .= 0
    end

    return solstep
end

function GeometricBase.reset!(solstep::SolutionStepDAE)
    for i in eachhistory(solstep)
        solstep.t̄[i]  = solstep.t̄[i-1]
        solstep.q̄[i] .= solstep.q̄[i-1]
        solstep.λ̄[i] .= solstep.λ̄[i-1]
        solstep.v̄[i] .= solstep.v̄[i-1]
        solstep.ū[i] .= solstep.ū[i-1]
    end

    return solstep
end

function update!(solstep::SolutionStepDAE{DT}, Δt, Δq::Vector{DT}, λ::Vector{DT}) where {DT}
    solstep.t += Δt
    for k in eachindex(Δq)
        solstep.q[k], solstep.q̃[k] = compensated_summation(Δq[k], solstep.q̄[1][k], solstep.q̃[k])
    end
    solstep.λ .= λ
    return solstep
end
