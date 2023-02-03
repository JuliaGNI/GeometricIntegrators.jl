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

    function SolutionStepPODE(t::TT, q::AT, p::AT; history = 2, internal::IT = NamedTuple()) where {DT, TT, AT <: AbstractArray{DT}, IT}
        @assert history ≥ 1

        t̄ = OffsetVector([zero(t) for _ in 1:history+1], 0:history)
        q̄ = OffsetVector([zero(q) for _ in 1:history+1], 0:history)
        p̄ = OffsetVector([zero(p) for _ in 1:history+1], 0:history)
        v̄ = OffsetVector([vectorfield(q) for _ in 1:history+1], 0:history)
        f̄ = OffsetVector([vectorfield(p) for _ in 1:history+1], 0:history)

        new{DT,TT,AT,typeof(v̄[0]),typeof(f̄[0]),IT,history}(t̄[0], q̄[0], p̄[0], v̄[0], f̄[0], t̄, q̄, p̄, v̄, f̄, zero(q), zero(p), internal)
    end
end

nhistory(::SolutionStepPODE{DT,TT,AT,VT,FT,IT,NT}) where {DT,TT,AT,VT,FT,IT,NT} = NT

current(solstep::SolutionStepPODE) = (t = solstep.t, q = solstep.q, p = solstep.p)
previous(solstep::SolutionStepPODE) = (t = solstep.t̄[1], q = solstep.q̄[1], p = solstep.p̄[1])
history(solstep::SolutionStepPODE) = (t = solstep.t̄, q = solstep.q̄, p = solstep.p̄)


function update_vector_fields!(solstep::SolutionStepPODE, problem::Union{PODEProblem,HODEProblem}, i=0)
    functions(problem).v(solstep.v̄[i], solstep.t̄[i], solstep.q̄[i], solstep.p̄[i])
    functions(problem).f(solstep.f̄[i], solstep.t̄[i], solstep.q̄[i], solstep.p̄[i])
end

function update_vector_fields!(solstep::SolutionStepPODE, problem::Union{IODEProblem,LODEProblem}, i=0)
    functions(problem).v̄(solstep.v̄[i], solstep.t̄[i], solstep.q̄[i])
    functions(problem).f̄(solstep.f̄[i], solstep.t̄[i], solstep.q̄[i], solstep.v̄[i])
    functions(problem).ϑ(solstep.p̄[i], solstep.t̄[i], solstep.q̄[i], solstep.v̄[i])
end

function initialize!(solstep::SolutionStepPODE, problem::AbstractProblemPODE, extrap::Extrapolation = default_extrapolation())
    solstep.t  = initial_conditions(problem).t
    solstep.q .= initial_conditions(problem).q
    solstep.p .= initial_conditions(problem).p

    solstep.q̃ .= 0
    solstep.p̃ .= 0

    update_vector_fields!(solstep, problem)

    for i in eachhistory(solstep)
        solstep.t̄[i] = solstep.t̄[i-1] - timestep(problem)
        extrapolate!(solstep.t̄[i-1], solstep.q̄[i-1], solstep.p̄[i-1], solstep.t̄[i], solstep.q̄[i], solstep.p̄[i], problem, extrap)
        update_vector_fields!(solstep, problem, i)
    end

    return solstep
end

function update!(solstep::SolutionStepPODE{DT,TT,AT}, Δq::AT, Δp::AT) where {DT,TT,AT}
    for k in eachindex(Δq,Δp)
        solstep.q[k], solstep.q̃[k] = compensated_summation(Δq[k], solstep.q̄[1][k], solstep.q̃[k])
        solstep.p[k], solstep.p̃[k] = compensated_summation(Δp[k], solstep.p̄[1][k], solstep.p̃[k])
    end
    return solstep
end

function GeometricBase.reset!(solstep::SolutionStepPODE, Δt)
    for i in backwardhistory(solstep)
        solstep.t̄[i]  = solstep.t̄[i-1]
        solstep.q̄[i] .= solstep.q̄[i-1]
        solstep.p̄[i] .= solstep.p̄[i-1]
        solstep.v̄[i] .= solstep.v̄[i-1]
        solstep.f̄[i] .= solstep.f̄[i-1]
    end

    solstep.t = solstep.t̄[0] += Δt

    return solstep
end


function Base.copy!(solstep::SolutionStepPODE, sol::NamedTuple)
    solstep.t  = sol.t
    solstep.q .= sol.q
    solstep.p .= sol.p
    solstep.q̃ .= 0
    solstep.p̃ .= 0
    solstep.v .= 0
    solstep.f .= 0

    # for i in backwardhistory(solstep)
    #     solstep.t̄[i]  = 0
    #     solstep.q̄[i] .= 0
    #     solstep.p̄[i] .= 0
    #     solstep.v̄[i] .= 0
    #     solstep.f̄[i] .= 0
    # end

    return solstep
end
