"""
Solution step for an [`ODEProblem`](@ref).

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
* `v`: vector field of q
* `v̄`: vector fields of q̄
* `internal`: internal variables of the integrator (e.g., internal stages of a Runge-Kutta methods or solver output)

### Constructors

```julia
SolutionStepODE(t::TT, q::AT, internal::IT=NamedTuple())
```

"""
mutable struct SolutionStepODE{
            DT <: Number, 
            TT <: Real, 
            AT <: AbstractArray{DT}, 
            VT <: AbstractArray{DT}, 
            IT <: NamedTuple,
            NT} <: SolutionStep{DT,TT,AT}

    t::TT
    q::AT
    v::VT

    t̄::OffsetVector{TT, Vector{TT}}
    q̄::OffsetVector{AT, Vector{AT}}
    v̄::OffsetVector{VT, Vector{VT}}

    q̃::AT

    internal::IT

    function SolutionStepODE(t::TT, q::AT; history = 2, internal::IT = NamedTuple()) where {DT, TT, AT <: AbstractArray{DT}, IT}
        @assert history ≥ 1

        t̄ = OffsetVector([zero(t) for _ in 1:history+1], 0:history)
        q̄ = OffsetVector([zero(q) for _ in 1:history+1], 0:history)
        v̄ = OffsetVector([vectorfield(q) for _ in 1:history+1], 0:history)

        new{DT,TT,AT,typeof(v̄[0]),IT,history}(t̄[0], q̄[0], v̄[0], t̄, q̄, v̄, zero(q), internal)
    end
end


nhistory(::SolutionStepODE{DT,TT,AT,VT,IT,NT}) where {DT,TT,AT,VT,IT,NT} = NT

current(solstep::SolutionStepODE) = (t = solstep.t, q = solstep.q)
previous(solstep::SolutionStepODE) = (t = solstep.t̄[1], q = solstep.q̄[1])
history(solstep::SolutionStepODE) = (t = solstep.t̄, q = solstep.q̄)


function update_vector_fields!(solstep::SolutionStepODE, problem::AbstractProblemODE, i=0)
    functions(problem).v(solstep.v̄[i], solstep.t̄[i], solstep.q̄[i])
end

# function update_vector_fields!(solstep::SolutionStepODE, problem::SODEProblem, i=0)
#     # TODO: add proper implementation !!!
# end

function initialize!(solstep::SolutionStepODE, problem::AbstractProblemODE, extrap::Extrapolation = default_extrapolation())
    solstep.t  = initial_conditions(problem).t
    solstep.q .= initial_conditions(problem).q
    solstep.q̃ .= 0

    update_vector_fields!(solstep, problem)

    for i in eachhistory(solstep)
        solstep.t̄[i] = solstep.t̄[i-1] - timestep(problem)
        extrapolate!(solstep.t̄[i-1], solstep.q̄[i-1], solstep.t̄[i], solstep.q̄[i], problem, extrap)
        update_vector_fields!(solstep, problem, i)
    end

    return solstep
end

function initialize!(solstep::SolutionStepODE, problem::SODEProblem, extrap::Extrapolation = default_extrapolation())
    solstep.t  = initial_conditions(problem).t
    solstep.q .= initial_conditions(problem).q
    solstep.q̃ .= 0

    return solstep
end

function update!(solstep::SolutionStepODE, Δq)
    for k in eachindex(Δq)
        solstep.q[k], solstep.q̃[k] = compensated_summation(solstep.q̄[1][k], Δq[k], solstep.q̃[k])
    end
    return solstep
end

function GeometricBase.reset!(solstep::SolutionStepODE, Δt)
    for i in backwardhistory(solstep)
        solstep.t̄[i]  = solstep.t̄[i-1]
        solstep.q̄[i] .= solstep.q̄[i-1]
        solstep.v̄[i] .= solstep.v̄[i-1]
    end

    solstep.t = solstep.t̄[0] += Δt

    return solstep
end


function Base.copy!(solstep::SolutionStepODE, sol::NamedTuple)
    solstep.t  = sol.t
    solstep.q .= sol.q
    solstep.q̃ .= 0
    solstep.v .= 0

    # for i in backwardhistory(solstep)
    #     solstep.t̄[i]  = 0
    #     solstep.q̄[i] .= 0
    #     solstep.v̄[i] .= 0
    # end

    return solstep
end
