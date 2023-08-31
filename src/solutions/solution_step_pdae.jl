"""
Solution step for a [`PDAEProblem`](@ref).

### Parameters

* `DT`: data type
* `TT`: time step type
* `AT`: array type
* `IT`: internal variable types

### Fields

* `t`: time of current time step
* `t̄`: time of previous time steps
* `q`: current solution of q
* `p`: current solution of p
* `λ`: current solution of λ
* `q̄`: previous solution of q
* `p̄`: previous solution of p
* `λ̄`: previous solution of λ
* `v`: vector field of q
* `f`: vector field of p
* `v̄`: vector field of q̄
* `f̄`: vector field of p̄
* `u`: projective vector field of q
* `ū`: projective vector field of q̄
* `g`: projective vector field of p
* `ḡ`: projective vector field of p̄
* `q̃`: compensated summation error of q
* `p̃`: compensated summation error of p
* `internal`: internal variables of the integrator (e.g., internal stages of a Runge-Kutta methods or solver output)

### Constructors

```julia
SolutionStepPDAE(t::TT, q::AT, p::AT, λ::AT; nhistory=1, internal::IT=NamedTuple())
```

"""
mutable struct SolutionStepPDAE{
                DT <: Number,
                TT <: Real,
                AT <: AbstractArray{DT},
                ΛT <: AbstractArray{DT},
                VT <: AbstractArray{DT},
                FT <: AbstractArray{DT},
                HT <: NamedTuple,
                IT <: NamedTuple,
                NT} <: SolutionStep{DT,TT,AT}
    
    q::AT
    p::AT
    λ::ΛT
    v::VT
    f::FT
    u::VT
    g::FT

    q̄::AT
    p̄::AT
    λ̄::ΛT
    v̄::VT
    f̄::FT
    ū::VT
    ḡ::FT

    q̃::AT
    p̃::AT

    history::HT
    internal::IT

    function SolutionStepPDAE(t::TT, q::AT, p::AT, λ::ΛT; nhistory = 2, internal::IT = NamedTuple()) where {DT, TT, AT <: AbstractArray{DT}, ΛT <: AbstractArray{DT}, IT}
        # TODO: nhistory should default to 1 and set to higher values by integrator / initial guess method
        @assert nhistory ≥ 1

        history = (
            t = OffsetVector([zero(t) for _ in 0:nhistory], 0:nhistory),
            q = OffsetVector([zero(q) for _ in 0:nhistory], 0:nhistory),
            p = OffsetVector([zero(p) for _ in 0:nhistory], 0:nhistory),
            λ = OffsetVector([zero(λ) for _ in 0:nhistory], 0:nhistory),
            v = OffsetVector([vectorfield(q) for _ in 0:nhistory], 0:nhistory),
            f = OffsetVector([vectorfield(p) for _ in 0:nhistory], 0:nhistory),
            u = OffsetVector([vectorfield(q) for _ in 0:nhistory], 0:nhistory),
            g = OffsetVector([vectorfield(p) for _ in 0:nhistory], 0:nhistory),
        )

        q = history.q[0]
        p = history.p[0]
        λ = history.λ[0]
        v = history.v[0]
        f = history.f[0]
        u = history.u[0]
        g = history.g[0]

        q̄ = history.q[1]
        p̄ = history.p[1]
        λ̄ = history.λ[1]
        v̄ = history.v[1]
        f̄ = history.f[1]
        ū = history.u[1]
        ḡ = history.g[1]

        q̃ = zero(q)
        p̃ = zero(p)

        new{DT, TT, AT, ΛT, typeof(v), typeof(f), typeof(history), IT, nhistory}(q, p, λ, v, f, u, g, q̄, p̄, λ̄, v̄, f̄, ū, ḡ, q̃, p̃, history, internal)
    end
end

@inline function Base.hasproperty(::SolutionStepPDAE, s::Symbol)
    s == :t || s == :t̄ || hasfield(SolutionStepPDAE, s)
end

@inline function Base.getproperty(solstep::SolutionStepPDAE, s::Symbol)
    if s == :t
        return history(solstep).t[0]
    elseif s == :t̄
        return history(solstep).t[1]
    else
        return getfield(solstep, s)
    end
end

@inline function Base.setproperty!(solstep::SolutionStepPDAE, s::Symbol, val)
    if s == :t
        return history(solstep).t[0] = val
    elseif s == :t̄
        return history(solstep).t[1] = val
    else
        return setfield!(solstep, s, val)
    end
end

nhistory(::SolutionStepPDAE{DT,TT,AT,ΛT,VT,FT,HT,IT,NT}) where {DT,TT,AT,ΛT,VT,FT,HT,IT,NT} = NT

current(solstep::SolutionStepPDAE) = (t = solstep.t, q = solstep.q, p = solstep.p, λ = solstep.λ)
previous(solstep::SolutionStepPDAE) = (t = solstep.t̄, q = solstep.q̄, p = solstep.p̄, λ = solstep.λ̄)
history(solstep::SolutionStepPDAE) = solstep.history
history(solstep::SolutionStepPDAE, i::Int) = (
    t = history(solstep).t[i],
    q = history(solstep).q[i],
    p = history(solstep).p[i],
    λ = history(solstep).λ[i],
    v = history(solstep).v[i],
    f = history(solstep).f[i],
    u = history(solstep).u[i],
    g = history(solstep).g[i])


function update_vector_fields!(solstep::SolutionStepPDAE, problem::Union{PDAEProblem,HDAEProblem}, i=0)
    functions(problem).v(history(solstep).v[i], history(solstep).t[i], history(solstep).q[i], history(solstep).p[i])
    functions(problem).f(history(solstep).f[i], history(solstep).t[i], history(solstep).q[i], history(solstep).p[i])
    functions(problem).u(history(solstep).u[i], history(solstep).t[i], history(solstep).q[i], history(solstep).p[i], history(solstep).λ[i])
    functions(problem).g(history(solstep).g[i], history(solstep).t[i], history(solstep).q[i], history(solstep).p[i], history(solstep).λ[i])
end

function update_vector_fields!(solstep::SolutionStepPDAE, problem::Union{IDAEProblem,LDAEProblem}, i=0)
    functions(problem).v̄(history(solstep).v[i], history(solstep).t[i], history(solstep).q[i], history(solstep).p[i])
    functions(problem).f̄(history(solstep).f[i], history(solstep).t[i], history(solstep).q[i], history(solstep).v[i])
    functions(problem).u(history(solstep).u[i], history(solstep).t[i], history(solstep).q[i], history(solstep).v[i], history(solstep).p[i], history(solstep).λ[i])
    functions(problem).g(history(solstep).g[i], history(solstep).t[i], history(solstep).q[i], history(solstep).v[i], history(solstep).p[i], history(solstep).λ[i])
end

update_implicit_functions!(::SolutionStepPDAE, ::Union{PDAEProblem,HDAEProblem}, i=0) = nothing

function update_implicit_functions!(solstep::SolutionStepPDAE, problem::Union{IDAEProblem,LDAEProblem}, i=0)
    functions(problem).ϑ(history(solstep).p[i], history(solstep).t[i], history(solstep).q[i], history(solstep).v[i])
end

function initialize!(solstep::SolutionStepPDAE, problem::Union{PDAEProblem, HDAEProblem, IDAEProblem, LDAEProblem}, extrap::Extrapolation = default_extrapolation())
    solstep.t  = initial_conditions(problem).t
    solstep.q .= initial_conditions(problem).q
    solstep.p .= initial_conditions(problem).p
    solstep.λ .= initial_conditions(problem).λ

    solstep.q̃ .= 0
    solstep.p̃ .= 0

    update_vector_fields!(solstep, problem)
    # update_implicit_functions!(solstep, problem)

    for i in eachhistory(solstep)
        history(solstep).t[i] = history(solstep).t[i-1] - timestep(problem)
        extrapolate!(history(solstep).t[i-1], history(solstep).q[i-1], history(solstep).p[i-1], history(solstep).t[i], history(solstep).q[i], history(solstep).p[i], problem, extrap)
        update_vector_fields!(solstep, problem, i)
        update_implicit_functions!(solstep, problem, i)
    end

    return solstep
end

function update!(solstep::SolutionStepPDAE{DT,TT,AT,ΛT}, Δq::AT, Δp::AT) where {DT,TT,AT,ΛT}
    for k in eachindex(Δq,Δp)
        solstep.q[k], solstep.q̃[k] = compensated_summation(Δq[k], solstep.q[k], solstep.q̃[k])
        solstep.p[k], solstep.p̃[k] = compensated_summation(Δp[k], solstep.p[k], solstep.p̃[k])
    end
    return solstep
end

function update!(solstep::SolutionStepPDAE{DT,TT,AT,ΛT}, q̇::AT, ṗ::AT, Δt::TT) where {DT,TT,AT,ΛT}
    for k in eachindex(q̇,ṗ)
        solstep.q[k], solstep.q̃[k] = compensated_summation(Δt * q̇[k], solstep.q[k], solstep.q̃[k])
        solstep.p[k], solstep.p̃[k] = compensated_summation(Δt * ṗ[k], solstep.p[k], solstep.p̃[k])
    end
    return solstep
end

function update!(solstep::SolutionStepPDAE{DT,TT,AT,ΛT}, Δq::AT, Δp::AT, λ::ΛT) where {DT,TT,AT,ΛT}
    update!(solstep, Δq, Δp)
    copyto!(solstep.λ, λ)
    return solstep
end

function update!(solstep::SolutionStepPDAE{DT,TT,AT,ΛT}, q̇::AT, ṗ::AT, λ::ΛT, Δt::TT) where {DT,TT,AT,ΛT}
    update!(solstep, q̇, ṗ, Δt)
    copyto!(solstep.λ, λ)
    return solstep
end

function GeometricBase.reset!(solstep::SolutionStepPDAE, Δt)
    for i in backwardhistory(solstep)
        history(solstep).t[i]  = history(solstep).t[i-1]
        history(solstep).q[i] .= history(solstep).q[i-1]
        history(solstep).p[i] .= history(solstep).p[i-1]
        history(solstep).λ[i] .= history(solstep).λ[i-1]
        history(solstep).v[i] .= history(solstep).v[i-1]
        history(solstep).f[i] .= history(solstep).f[i-1]
        history(solstep).u[i] .= history(solstep).u[i-1]
        history(solstep).g[i] .= history(solstep).g[i-1]
    end

    solstep.t = solstep.t̄ + Δt

    return solstep
end


function Base.copy!(solstep::SolutionStepPDAE, sol::NamedTuple)
    solstep.t  = sol.t
    solstep.q .= sol.q
    solstep.p .= sol.p
    solstep.λ .= sol.λ
    solstep.q̃ .= 0
    solstep.p̃ .= 0
    solstep.v .= 0
    solstep.f .= 0
    solstep.u .= 0
    solstep.g .= 0

    # for i in backwardhistory(solstep)
    #     solstep.t̄[i]  = 0
    #     solstep.q̄[i] .= 0
    #     solstep.p̄[i] .= 0
    #     solstep.λ̄[i] .= 0
    #     solstep.v̄[i] .= 0
    #     solstep.f̄[i] .= 0
    #     solstep.ū[i] .= 0
    #     solstep.ḡ[i] .= 0
    # end

    return solstep
end
