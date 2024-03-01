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
* `λ`: current solution of λ
* `μ`: current solution of μ
* `q̄`: previous solutions of q
* `λ̄`: previous solutions of λ
* `μ̄`: previous solution of μ
* `v`: vector field of q
* `v̄`: vector field of q̄
* `u`: projective vector field of q
* `ū`: projective vector field of q̄
* `q̃`: compensated summation error of q
* `internal`: internal variables of the integrator (e.g., internal stages of a Runge-Kutta methods or solver output)

### Constructors

```julia
SolutionStepDAE(t::TT, q::AT, λ::AT, internal::IT=NamedTuple())
```

"""
struct SolutionStepDAE{
            DT <: Number,
            TT <: Real,
            AT <: AbstractArray{DT},
            ΛT <: AbstractArray{DT},
            VT <: AbstractArray{DT},
            HT <: NamedTuple,
            IT <: NamedTuple,
            paramsType <: OptionalParameters,
            NT} <: SolutionStep{DT,TT,AT}

    q::AT
    λ::ΛT
    μ::ΛT
    v::VT
    u::VT

    q̄::AT
    λ̄::ΛT
    μ̄::ΛT
    v̄::VT
    ū::VT

    q̃::AT

    history::HT
    internal::IT

    parameters::paramsType

    function SolutionStepDAE(t::TT, q::AT, λ::ΛT, μ::ΛT, parameters; nhistory = 2, internal::IT = NamedTuple()) where {DT, TT, AT <: AbstractArray{DT}, ΛT <: AbstractArray{DT}, IT}
        # TODO: nhistory should default to 1 and set to higher values by integrator / initial guess method
        @assert nhistory ≥ 1

        history = (
            t = OffsetVector([zero(t) for _ in 0:nhistory], 0:nhistory),
            q = OffsetVector([zero(q) for _ in 0:nhistory], 0:nhistory),
            λ = OffsetVector([zero(λ) for _ in 0:nhistory], 0:nhistory),
            μ = OffsetVector([zero(μ) for _ in 0:nhistory], 0:nhistory),
            v = OffsetVector([vectorfield(q) for _ in 0:nhistory], 0:nhistory),
            u = OffsetVector([vectorfield(q) for _ in 0:nhistory], 0:nhistory),
        )

        q = history.q[0]
        λ = history.λ[0]
        μ = history.μ[0]
        v = history.v[0]
        u = history.u[0]

        q̄ = history.q[1]
        λ̄ = history.λ[1]
        μ̄ = history.μ[1]
        v̄ = history.v[1]
        ū = history.u[1]

        q̃ = zero(q)

        new{DT, TT, AT, ΛT, typeof(v), typeof(history), IT, typeof(parameters), nhistory}(q, λ, μ, v, u, q̄, λ̄, μ̄, v̄, ū, q̃, history, internal, parameters)
    end
end

@inline function Base.hasproperty(::SolutionStepDAE, s::Symbol)
    s == :t || s == :t̄ || hasfield(SolutionStepDAE, s)
end

@inline function Base.getproperty(solstep::SolutionStepDAE, s::Symbol)
    if s == :t
        return history(solstep).t[0]
    elseif s == :t̄
        return history(solstep).t[1]
    else
        return getfield(solstep, s)
    end
end

@inline function Base.setproperty!(solstep::SolutionStepDAE, s::Symbol, val)
    if s == :t
        return history(solstep).t[0] = val
    elseif s == :t̄
        return history(solstep).t[1] = val
    else
        return setfield!(solstep, s, val)
    end
end

nhistory(::SolutionStepDAE{DT,TT,AT,ΛT,VT,HT,IT,PT,NT}) where {DT,TT,AT,ΛT,VT,HT,IT,PT,NT} = NT

current(solstep::SolutionStepDAE) = (t = solstep.t, q = solstep.q, λ = solstep.λ, μ = solstep.μ)
previous(solstep::SolutionStepDAE) = (t = solstep.t̄, q = solstep.q̄, λ = solstep.λ̄, μ = solstep.μ̄)
history(solstep::SolutionStepDAE) = solstep.history
history(solstep::SolutionStepDAE, i::Int) = (
    t = history(solstep).t[i],
    q = history(solstep).q[i],
    λ = history(solstep).λ[i],
    μ = history(solstep).μ[i],
    v = history(solstep).v[i],
    u = history(solstep).u[i])
internal(solstep::SolutionStepDAE) = solstep.internal
parameters(solstep::SolutionStepDAE) = solstep.parameters


function update_vector_fields!(solstep::SolutionStepDAE, problem::DAEProblem, i=0)
    functions(problem).v(history(solstep).v[i], history(solstep).t[i], history(solstep).q[i], parameters(problem))
    functions(problem).u(history(solstep).u[i], history(solstep).t[i], history(solstep).q[i], history(solstep).λ[i], parameters(problem))
end

function initialize!(solstep::SolutionStepDAE, problem::DAEProblem, extrap::Extrapolation = default_extrapolation())
    solstep.t  = initial_conditions(problem).t
    solstep.q .= initial_conditions(problem).q
    solstep.λ .= initial_conditions(problem).λ
    solstep.μ .= initial_conditions(problem).μ
    solstep.q̃ .= 0

    update_vector_fields!(solstep, problem)

    for i in eachhistory(solstep)
        history(solstep).t[i] = history(solstep).t[i-1] - timestep(problem)
        extrapolate!(history(solstep).t[i-1], history(solstep).q[i-1], history(solstep).t[i], history(solstep).q[i], problem, extrap)
        update_vector_fields!(solstep, problem, i)
    end

    return solstep
end

function update!(solstep::SolutionStepDAE{DT,TT,AT,ΛT}, Δq::Union{AT,AbstractArray}) where {DT,TT,AT,ΛT}
    for k in eachindex(Δq)
        solstep.q[k], solstep.q̃[k] = compensated_summation(Δq[k], solstep.q[k], solstep.q̃[k])
    end
    return solstep
end

function update!(solstep::SolutionStepDAE{DT,TT,AT,ΛT}, q̇::Union{AT,AbstractArray}, Δt::TT) where {DT,TT,AT,ΛT}
    for k in eachindex(q̇)
        solstep.q[k], solstep.q̃[k] = compensated_summation(Δt * q̇[k], solstep.q[k], solstep.q̃[k])
    end
    return solstep
end

function update!(solstep::SolutionStepDAE{DT,TT,AT,ΛT}, Δq::Union{AT,AbstractArray}, λ::Union{ΛT,AbstractArray}) where {DT,TT,AT,ΛT}
    update!(solstep, Δq)
    copyto!(solstep.λ, λ)
    return solstep
end

function update!(solstep::SolutionStepDAE{DT,TT,AT,ΛT}, Δq::Union{AT,AbstractArray}, λ::Union{ΛT,AbstractArray}, μ::Union{ΛT,AbstractArray}) where {DT,TT,AT,ΛT}
    update!(solstep, Δq, λ)
    copyto!(solstep.μ, μ)
    return solstep
end

function update!(solstep::SolutionStepDAE{DT,TT,AT,ΛT}, q̇::Union{AT,AbstractArray}, λ::Union{ΛT,AbstractArray}, Δt::TT) where {DT,TT,AT,ΛT}
    update!(solstep, q̇, Δt)
    copyto!(solstep.λ, λ)
    return solstep
end

function update!(solstep::SolutionStepDAE{DT,TT,AT,ΛT}, q̇::Union{AT,AbstractArray}, λ::Union{ΛT,AbstractArray}, μ::Union{ΛT,AbstractArray}, Δt::TT) where {DT,TT,AT,ΛT}
    update!(solstep, q̇, λ, Δt)
    copyto!(solstep.μ, μ)
    return solstep
end

function GeometricBase.reset!(solstep::SolutionStepDAE, Δt)
    for i in backwardhistory(solstep)
        history(solstep).t[i]  = history(solstep).t[i-1]
        history(solstep).q[i] .= history(solstep).q[i-1]
        history(solstep).λ[i] .= history(solstep).λ[i-1]
        history(solstep).μ[i] .= history(solstep).μ[i-1]
        history(solstep).v[i] .= history(solstep).v[i-1]
        history(solstep).u[i] .= history(solstep).u[i-1]
    end

    solstep.t = solstep.t̄ + Δt

    return solstep
end


function Base.copy!(solstep::SolutionStepDAE, sol::NamedTuple)
    solstep.t  = sol.t
    solstep.q .= sol.q
    solstep.λ .= sol.λ
    solstep.μ .= sol.μ
    solstep.q̃ .= 0
    solstep.v .= 0

    # for i in backwardhistory(solstep)
    #     solstep.t̄[i]  = 0
    #     solstep.q̄[i] .= 0
    #     solstep.λ̄[i] .= 0
    #     solstep.μ̄[i] .= 0
    #     solstep.v̄[i] .= 0
    #     solstep.ū[i] .= 0
    # end

    return solstep
end
