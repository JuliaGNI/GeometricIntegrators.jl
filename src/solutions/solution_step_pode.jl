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
* `p`: current solution of p
* `q̄`: previous solution of q
* `p̄`: previous solution of p
* `v`: vector field of q
* `f`: vector field of p
* `v̄`: vector field of q̄
* `f̄`: vector field of p̄
* `q̃`: compensated summation error of q
* `p̃`: compensated summation error of p
* `internal`: internal variables of the integrator (e.g., internal stages of a Runge-Kutta methods or solver output)

### Constructors

```julia
SolutionStepPODE(t::TT, q::AT, p::AT; nhistory=1, internal::IT=NamedTuple())
```

"""
mutable struct SolutionStepPODE{
            DT <: Number,
            TT <: Real,
            AT <: AbstractArray{DT},
            VT <: AbstractArray{DT},
            FT <: AbstractArray{DT},
            HT <: NamedTuple,
            IT <: NamedTuple,
            NT} <: SolutionStep{DT,TT,AT}

    q::AT
    p::AT
    v::VT
    f::FT

    q̄::AT
    p̄::AT
    v̄::VT
    f̄::FT

    q̃::AT
    p̃::AT

    history::HT
    internal::IT

    function SolutionStepPODE(t::TT, q::AT, p::AT; nhistory = 2, internal::IT = NamedTuple()) where {DT, TT, AT <: AbstractArray{DT}, IT}
        # TODO: nhistory should default to 1 and set to higher values by integrator / initial guess method
        @assert nhistory ≥ 1

        history = (
            t = OffsetVector([zero(t) for _ in 0:nhistory], 0:nhistory),
            q = OffsetVector([zero(q) for _ in 0:nhistory], 0:nhistory),
            p = OffsetVector([zero(p) for _ in 0:nhistory], 0:nhistory),
            v = OffsetVector([vectorfield(q) for _ in 0:nhistory], 0:nhistory),
            f = OffsetVector([vectorfield(p) for _ in 0:nhistory], 0:nhistory),
        )

        q = history.q[0]
        p = history.p[0]
        v = history.v[0]
        f = history.f[0]

        q̄ = history.q[1]
        p̄ = history.p[1]
        v̄ = history.v[1]
        f̄ = history.f[1]

        q̃ = zero(q)
        p̃ = zero(p)

        new{DT, TT, AT, typeof(v), typeof(f), typeof(history), IT, nhistory}(q, p, v, f, q̄, p̄, v̄, f̄, q̃, p̃, history, internal)
    end
end

@inline function Base.hasproperty(::SolutionStepPODE, s::Symbol)
    s == :t || s == :t̄ || hasfield(SolutionStepPODE, s)
end

@inline function Base.getproperty(solstep::SolutionStepPODE, s::Symbol)
    if s == :t
        return history(solstep).t[0]
    elseif s == :t̄
        return history(solstep).t[1]
    else
        return getfield(solstep, s)
    end
end

@inline function Base.setproperty!(solstep::SolutionStepPODE, s::Symbol, val)
    if s == :t
        return history(solstep).t[0] = val
    elseif s == :t̄
        return history(solstep).t[1] = val
    else
        return setfield!(solstep, s, val)
    end
end

nhistory(::SolutionStepPODE{DT,TT,AT,VT,FT,HT,IT,NT}) where {DT,TT,AT,VT,FT,HT,IT,NT} = NT

current(solstep::SolutionStepPODE) = (t = solstep.t, q = solstep.q, p = solstep.p)
previous(solstep::SolutionStepPODE) = (t = solstep.t̄, q = solstep.q̄, p = solstep.p̄)
history(solstep::SolutionStepPODE) = solstep.history
history(solstep::SolutionStepPODE, i::Int) = (
    t = history(solstep).t[i],
    q = history(solstep).q[i],
    p = history(solstep).p[i],
    v = history(solstep).v[i],
    f = history(solstep).f[i])


function update_vector_fields!(solstep::SolutionStepPODE, problem::Union{PODEProblem,HODEProblem}, i=0)
    functions(problem).v(history(solstep).v[i], history(solstep).t[i], history(solstep).q[i], history(solstep).p[i])
    functions(problem).f(history(solstep).f[i], history(solstep).t[i], history(solstep).q[i], history(solstep).p[i])
end

function update_vector_fields!(solstep::SolutionStepPODE, problem::Union{IODEProblem,LODEProblem}, i=0)
    functions(problem).v̄(history(solstep).v[i], history(solstep).t[i], history(solstep).q[i])
    functions(problem).f̄(history(solstep).f[i], history(solstep).t[i], history(solstep).q[i], history(solstep).v[i])
end

update_implicit_functions!(::SolutionStepPODE, ::Union{PODEProblem,HODEProblem}, i=0) = nothing

function update_implicit_functions!(solstep::SolutionStepPODE, problem::Union{IODEProblem,LODEProblem}, i=0)
    functions(problem).ϑ(history(solstep).p[i], history(solstep).t[i], history(solstep).q[i], history(solstep).v[i])
end

function initialize!(solstep::SolutionStepPODE, problem::AbstractProblemPODE, extrap::Extrapolation = default_extrapolation())
    solstep.t  = initial_conditions(problem).t
    solstep.q .= initial_conditions(problem).q
    solstep.p .= initial_conditions(problem).p

    solstep.q̃ .= 0
    solstep.p̃ .= 0

    update_vector_fields!(solstep, problem)

    for i in eachhistory(solstep)
        history(solstep).t[i] = history(solstep).t[i-1] - timestep(problem)
        extrapolate!(history(solstep).t[i-1], history(solstep).q[i-1], history(solstep).p[i-1], history(solstep).t[i], history(solstep).q[i], history(solstep).p[i], problem, extrap)
        update_vector_fields!(solstep, problem, i)
        update_implicit_functions!(solstep, problem, i)
    end

    return solstep
end

function update!(solstep::SolutionStepPODE, Δq, Δp)
    for k in eachindex(Δq,Δp)
        solstep.q[k], solstep.q̃[k] = compensated_summation(Δq[k], solstep.q[k], solstep.q̃[k])
        solstep.p[k], solstep.p̃[k] = compensated_summation(Δp[k], solstep.p[k], solstep.p̃[k])
    end
    return solstep
end

function update!(solstep::SolutionStepPODE, y, z, Δt)
    for k in eachindex(y,z)
        solstep.q[k], solstep.q̃[k] = compensated_summation(Δt * y[k], solstep.q[k], solstep.q̃[k])
        solstep.p[k], solstep.p̃[k] = compensated_summation(Δt * z[k], solstep.p[k], solstep.p̃[k])
    end
    return solstep
end

function GeometricBase.reset!(solstep::SolutionStepPODE, Δt)
    for i in backwardhistory(solstep)
        history(solstep).t[i]  = history(solstep).t[i-1]
        history(solstep).q[i] .= history(solstep).q[i-1]
        history(solstep).p[i] .= history(solstep).p[i-1]
        history(solstep).v[i] .= history(solstep).v[i-1]
        history(solstep).f[i] .= history(solstep).f[i-1]
    end

    solstep.t = solstep.t̄ + Δt

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