"""
Solution step for an [`ODEProblem`](@ref).

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
* `v`: vector field of q
* `v̄`: vector field of q̄
* `q̃`: compensated summation error of q
* `internal`: internal variables of the integrator (e.g., internal stages of a Runge-Kutta methods or solver output)

### Constructors

```julia
SolutionStepODE(t::TT, q::AT; nhistory=1, internal::IT=NamedTuple())
```

"""
struct SolutionStepODE{
            DT <: Number,
            TT <: Real,
            AT <: AbstractArray{DT},
            VT <: AbstractArray{DT},
            HT <: NamedTuple,
            IT <: NamedTuple,
            paramsType <: OptionalParameters,
            NT} <: SolutionStep{DT,TT,AT}

    q::AT
    v::VT

    q̄::AT
    v̄::VT

    q̃::AT

    history::HT
    internal::IT

    parameters::paramsType

    function SolutionStepODE(t::TT, q::AT, parameters; nhistory = 2, internal::IT = NamedTuple()) where {DT, TT, AT <: AbstractArray{DT}, IT}
        # TODO: nhistory should default to 1 and set to higher values by integrator / initial guess method
        @assert nhistory ≥ 1

        history = (
            t = OffsetVector([zero(t) for _ in 0:nhistory], 0:nhistory),
            q = OffsetVector([zero(q) for _ in 0:nhistory], 0:nhistory),
            v = OffsetVector([vectorfield(q) for _ in 0:nhistory], 0:nhistory),
        )

        q = history.q[0]
        v = history.v[0]

        q̄ = history.q[1]
        v̄ = history.v[1]

        q̃ = zero(q)

        new{DT, TT, AT, typeof(v), typeof(history), IT, typeof(parameters), nhistory}(q, v, q̄, v̄, q̃, history, internal, parameters)
    end
end

@inline function Base.hasproperty(::SolutionStepODE, s::Symbol)
    s == :t || s == :t̄ || hasfield(SolutionStepODE, s)
end

@inline function Base.getproperty(solstep::SolutionStepODE, s::Symbol)
    if s == :t
        return history(solstep).t[0]
    elseif s == :t̄
        return history(solstep).t[1]
    else
        return getfield(solstep, s)
    end
end

@inline function Base.setproperty!(solstep::SolutionStepODE, s::Symbol, val)
    if s == :t
        return history(solstep).t[0] = val
    elseif s == :t̄
        return history(solstep).t[1] = val
    else
        return setfield!(solstep, s, val)
    end
end

nhistory(::SolutionStepODE{DT,TT,AT,VT,HT,IT,PT,NT}) where {DT,TT,AT,VT,HT,IT,PT,NT} = NT

current(solstep::SolutionStepODE) = (t = solstep.t, q = solstep.q)
previous(solstep::SolutionStepODE) = (t = solstep.t̄, q = solstep.q̄)
history(solstep::SolutionStepODE) = solstep.history
history(solstep::SolutionStepODE, i::Int) = (
    t = history(solstep).t[i],
    q = history(solstep).q[i],
    v = history(solstep).v[i])
internal(solstep::SolutionStepODE) = solstep.internal
parameters(solstep::SolutionStepODE) = solstep.parameters

function update_vector_fields!(solstep::SolutionStepODE, problem::Union{ODEProblem, SubstepProblem}, i=0)
    functions(problem).v(history(solstep).v[i], history(solstep).t[i], history(solstep).q[i], parameters(problem))
end

function update_vector_fields!(solstep::SolutionStepODE, problem::SODEProblem, i=0)
    initialguess(problem).v(history(solstep).v[i], history(solstep).t[i], history(solstep).q[i], parameters(problem))
end

function update_vector_fields!(solstep::SolutionStepODE, problem::DELEProblem, i=0)
    nothing
end

function initialize!(solstep::SolutionStepODE, sol::NamedTuple, problem::Union{ODEProblem, SODEProblem, SubstepProblem}, extrap::Extrapolation = default_extrapolation())
    solstep.t  = sol.t
    solstep.q .= sol.q
    solstep.q̃ .= 0

    update_vector_fields!(solstep, problem)

    for i in eachhistory(solstep)
        history(solstep).t[i] = history(solstep).t[i-1] - timestep(problem)
        soltmp = (t = history(solstep).t[i], q = history(solstep).q[i], v = history(solstep).v[i])
        hsttmp = (t = [history(solstep).t[i-1]], q = [history(solstep).q[i-1]], v = [history(solstep).v[i-1]])
        solutionstep!(soltmp, hsttmp, problem, extrap)
    end

    return solstep
end

function initialize!(solstep::SolutionStepODE, sol::NamedTuple, problem::DELEProblem, extrap::Extrapolation = default_extrapolation())
    solstep.t  = sol.t
    solstep.q .= sol.q

    history(solstep).t[1]  = sol.t - timestep(problem)
    history(solstep).q[1] .= sol.q̄

    return solstep
end

function initialize!(solstep::SolutionStepODE, problem::Union{ODEProblem, SODEProblem, SubstepProblem, DELEProblem}, args...)
    initialize!(solstep, initial_conditions(problem), problem, args...)
end

function update!(solstep::SolutionStepODE, Δq)
    for k in eachindex(Δq)
        solstep.q[k], solstep.q̃[k] = compensated_summation(Δq[k], solstep.q[k], solstep.q̃[k])
    end
    return solstep
end

function update!(solstep::SolutionStepODE, q̇, Δt)
    for k in eachindex(q̇)
        solstep.q[k], solstep.q̃[k] = compensated_summation(Δt * q̇[k], solstep.q[k], solstep.q̃[k])
    end
    return solstep
end

function GeometricBase.reset!(solstep::SolutionStepODE, Δt)
    for i in backwardhistory(solstep)
        history(solstep).t[i]  = history(solstep).t[i-1]
        history(solstep).q[i] .= history(solstep).q[i-1]
        history(solstep).v[i] .= history(solstep).v[i-1]
    end

    solstep.t = solstep.t̄ + Δt

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
