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
            HT <: NamedTuple,
            IT <: NamedTuple,
            NT} <: SolutionStep{DT,TT,AT}

    q::AT
    v::VT

    q̄::AT
    v̄::VT

    q̃::AT

    history::HT
    internal::IT

    function SolutionStepODE(t::TT, q::AT; nhistory = 2, internal::IT = NamedTuple()) where {DT, TT, AT <: AbstractArray{DT}, IT}
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

        new{DT, TT, AT, typeof(v), typeof(history), IT, nhistory}(q, v, q̄, v̄, q̃, history, internal)
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

nhistory(::SolutionStepODE{DT,TT,AT,VT,HT,IT,NT}) where {DT,TT,AT,VT,HT,IT,NT} = NT

current(solstep::SolutionStepODE) = (t = solstep.t, q = solstep.q)
previous(solstep::SolutionStepODE) = (t = solstep.t̄, q = solstep.q̄)
history(solstep::SolutionStepODE) = solstep.history
history(solstep::SolutionStepODE, i::Int) = (
    t = history(solstep).t[i],
    q = history(solstep).q[i],
    v = history(solstep).v[i])


function update_vector_fields!(solstep::SolutionStepODE, problem::AbstractProblemODE, i=0)
    functions(problem).v(history(solstep).v[i], history(solstep).t[i], history(solstep).q[i])
end

function update_vector_fields!(solstep::SolutionStepODE, problem::SODEProblem, i=0)
    # TODO: add proper implementation !!!
    # @warn "update_vector_fields!() method for SODEs is still missing"
end

function initialize!(solstep::SolutionStepODE, problem::AbstractProblemODE, extrap::Extrapolation = default_extrapolation())
    solstep.t  = initial_conditions(problem).t
    solstep.q .= initial_conditions(problem).q
    solstep.q̃ .= 0

    update_vector_fields!(solstep, problem)

    for i in eachhistory(solstep)
        history(solstep).t[i] = history(solstep).t[i-1] - timestep(problem)
        extrapolate!(history(solstep).t[i-1], history(solstep).q[i-1], history(solstep).t[i], history(solstep).q[i], problem, extrap)
        update_vector_fields!(solstep, problem, i)
    end

    return solstep
end

function initialize!(solstep::SolutionStepODE, problem::SODEProblem, extrap::Extrapolation = default_extrapolation())
    solstep.t  = initial_conditions(problem).t
    solstep.q .= initial_conditions(problem).q
    solstep.q̃ .= 0

    for i in eachhistory(solstep)
        history(solstep).t[i] = history(solstep).t[i-1] - timestep(problem)
        # TODO: add proper implementation !!!
        # extrapolate!(solstep.t̄[i-1], solstep.q̄[i-1], solstep.t̄[i], solstep.q̄[i], problem, extrap)
        # update_vector_fields!(solstep, problem, i)
    end

    return solstep
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
