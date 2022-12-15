"""
Solution step for a [`PODEProblem`](@ref).

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
* `p`: current solution of p
* `p̄`: previous solutions of p
* `p̃`: compensated summation error of p
* `v`: vector field of q
* `v̄`: vector fields of q̄
* `f`: vector field of p
* `f̄`: vector fields of p̄
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
            IT <: NamedTuple,
            NT} <: SolutionStep{DT,TT,AT}

    t::TT
    q::AT
    p::AT
    v::VT
    f::VT

    t̄::OffsetVector{TT, Vector{TT}}
    q̄::OffsetVector{AT, Vector{AT}}
    p̄::OffsetVector{AT, Vector{AT}}
    v̄::OffsetVector{VT, Vector{VT}}
    f̄::OffsetVector{VT, Vector{VT}}

    q̃::AT
    p̃::AT

    internal::IT

    function SolutionStepPODE(t::TT, q::AT, p::AT; history = 1, internal::IT = NamedTuple()) where {DT, TT, AT <: AbstractArray{DT}, IT}
        @assert history ≥ 1

        t̄ = OffsetVector([zero(t) for _ in 1:history+1], 0:history)
        q̄ = OffsetVector([zero(q) for _ in 1:history+1], 0:history)
        p̄ = OffsetVector([zero(p) for _ in 1:history+1], 0:history)
        v̄ = OffsetVector([vectorfield(q) for _ in 1:history+1], 0:history)
        f̄ = OffsetVector([vectorfield(p) for _ in 1:history+1], 0:history)

        t̄[0]  = t
        q̄[0] .= q
        p̄[0] .= p

        new{DT,TT,AT,typeof(v̄[0]),typeof(f̄[0]),IT,history}(t̄[0], q̄[0], p̄[0], v̄[0], f̄[0], t̄, q̄, p̄, v̄, f̄, zero(q), zero(p), internal)
    end
end

nhistory(::SolutionStepPODE{DT,TT,AT,VT,FT,IT,NT}) where {DT,TT,AT,VT,FT,IT,NT} = NT

current(solstep::SolutionStepPODE) = (t = solstep.t, q = solstep.q, p = solstep.p)
previous(solstep::SolutionStepPODE) = (t = solstep.t̄[1], q = solstep.q̄[1], p = solstep.p̄[1])
history(solstep::SolutionStepPODE) = (t = solstep.t̄, q = solstep.q̄, p = solstep.p̄)


function Base.copy!(solstep::SolutionStepPODE, sol::NamedTuple)
    solstep.t  = sol.t
    solstep.q .= sol.q
    solstep.p .= sol.p
    solstep.q̃ .= 0
    solstep.p̃ .= 0
    solstep.v .= 0
    solstep.f .= 0

    for i in eachhistory(solstep)
        solstep.t̄[i]  = 0
        solstep.q̄[i] .= 0
        solstep.p̄[i] .= 0
        solstep.v̄[i] .= 0
        solstep.f̄[i] .= 0
    end

    return solstep
end

function GeometricBase.reset!(solstep::SolutionStepPODE)
    for i in eachhistory(solstep)
        solstep.t̄[i]  = solstep.t̄[i-1]
        solstep.q̄[i] .= solstep.q̄[i-1]
        solstep.p̄[i] .= solstep.p̄[i-1]
        solstep.v̄[i] .= solstep.v̄[i-1]
        solstep.f̄[i] .= solstep.f̄[i-1]
    end

    return solstep
end

function update!(solstep::SolutionStepPODE{DT,TT,AT}, Δt::TT, Δq::AT, Δp::AT) where {DT,TT,AT}
    solstep.t += Δt
    for k in eachindex(Δq,Δp)
        solstep.q[k], solstep.q̃[k] = compensated_summation(Δq[k], solstep.q̄[1][k], solstep.q̃[k])
        solstep.p[k], solstep.p̃[k] = compensated_summation(Δp[k], solstep.p̄[1][k], solstep.p̃[k])
    end
    return solstep
end
