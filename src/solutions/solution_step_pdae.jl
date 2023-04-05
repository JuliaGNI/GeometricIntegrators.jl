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
* `q̄`: previous solutions of q
* `q̃`: compensated summation error of q
* `p`: current solution of p
* `p̄`: previous solutions of p
* `p̃`: compensated summation error of p
* `λ`: current solution of λ
* `λ̄`: previous solutions of λ
* `v`: vector field of q
* `v̄`: vector fields of q̄
* `f`: vector field of p
* `f̄`: vector fields of p̄
* `u`: projective vector field of q
* `ū`: projective vector fields of q̄
* `g`: projective vector field of p
* `ḡ`: projective vector fields of p̄
* `internal`: internal variables of the integrator (e.g., internal stages of a Runge-Kutta methods or solver output)

### Constructors

```julia
SolutionStepPDAE(t::TT, q::AT, p::AT, λ::AT, internal::IT=NamedTuple())
```

"""
mutable struct SolutionStepPDAE{
                DT <: Number,
                TT <: Real,
                AT <: AbstractArray{DT},
                ΛT <: AbstractArray{DT},
                VT <: AbstractArray{DT},
                FT <: AbstractArray{DT},
                IT <: NamedTuple,
                NT} <: SolutionStep{DT,TT,AT}
    
    t::TT
    q::AT
    p::AT
    λ::ΛT
    v::VT
    f::FT
    u::VT
    g::FT

    t̄::OffsetVector{TT, Vector{TT}}
    q̄::OffsetVector{AT, Vector{AT}}
    p̄::OffsetVector{AT, Vector{AT}}
    λ̄::OffsetVector{ΛT, Vector{ΛT}}
    v̄::OffsetVector{VT, Vector{VT}}
    f̄::OffsetVector{FT, Vector{FT}}
    ū::OffsetVector{VT, Vector{VT}}
    ḡ::OffsetVector{FT, Vector{FT}}

    q̃::AT
    p̃::AT

    internal::IT

    function SolutionStepPDAE(t::TT, q::AT, p::AT, λ::ΛT; history = 2, internal::IT = NamedTuple()) where {DT, TT, AT <: AbstractArray{DT}, ΛT <: AbstractArray{DT}, IT}
        @assert history ≥ 1

        t̄ = OffsetVector([zero(t) for _ in 1:history+1], 0:history)
        q̄ = OffsetVector([zero(q) for _ in 1:history+1], 0:history)
        p̄ = OffsetVector([zero(p) for _ in 1:history+1], 0:history)
        λ̄ = OffsetVector([zero(λ) for _ in 1:history+1], 0:history)
        v̄ = OffsetVector([vectorfield(q) for _ in 1:history+1], 0:history)
        f̄ = OffsetVector([vectorfield(p) for _ in 1:history+1], 0:history)
        ū = OffsetVector([vectorfield(q) for _ in 1:history+1], 0:history)
        ḡ = OffsetVector([vectorfield(p) for _ in 1:history+1], 0:history)

        new{DT,TT,AT,ΛT,typeof(v̄[0]),typeof(f̄[0]),IT,history}(t̄[0], q̄[0], p̄[0], λ̄[0], v̄[0], f̄[0], ū[0], ḡ[0], t̄, q̄, p̄, λ̄, v̄, f̄, ū, ḡ, zero(q), zero(p), internal)
    end
end


nhistory(::SolutionStepPDAE{DT,TT,AT,ΛT,VT,FT,IT,NT}) where {DT,TT,AT,ΛT,VT,FT,IT,NT} = NT

current(solstep::SolutionStepPDAE) = (t = solstep.t, q = solstep.q, p = solstep.p, λ = solstep.λ)
previous(solstep::SolutionStepPDAE) = (t = solstep.t̄[1], q = solstep.q̄[1], p = solstep.p̄[1], λ = solstep.λ̄[1])
history(solstep::SolutionStepPDAE) = (t = solstep.t̄, q = solstep.q̄, p = solstep.p̄, λ = solstep.λ̄)


function update_vector_fields!(solstep::SolutionStepPDAE, problem::Union{PDAEProblem,HDAEProblem}, i=0)
    functions(problem).v(solstep.v̄[i], solstep.t̄[i], solstep.q̄[i], solstep.p̄[i])
    functions(problem).f(solstep.f̄[i], solstep.t̄[i], solstep.q̄[i], solstep.p̄[i])
    functions(problem).u(solstep.ū[i], solstep.t̄[i], solstep.q̄[i], solstep.p̄[i], solstep.λ̄[i])
    functions(problem).g(solstep.ḡ[i], solstep.t̄[i], solstep.q̄[i], solstep.p̄[i], solstep.λ̄[i])
end

function update_vector_fields!(solstep::SolutionStepPDAE, problem::Union{IDAEProblem,LDAEProblem}, i=0)
    functions(problem).v̄(solstep.v̄[i], solstep.t̄[i], solstep.q̄[i])
    functions(problem).f̄(solstep.f̄[i], solstep.t̄[i], solstep.q̄[i], solstep.v̄[i])
    functions(problem).u(solstep.ū[i], solstep.t̄[i], solstep.q̄[i], solstep.v̄[i], solstep.p̄[i], solstep.λ̄[i])
    functions(problem).g(solstep.ḡ[i], solstep.t̄[i], solstep.q̄[i], solstep.v̄[i], solstep.p̄[i], solstep.λ̄[i])
end

update_implicit_functions!(::SolutionStepPDAE, ::Union{PDAEProblem,HDAEProblem}, i=0) = nothing

function update_implicit_functions!(solstep::SolutionStepPDAE, problem::Union{IDAEProblem,LDAEProblem}, i=0)
    functions(problem).ϑ(solstep.p̄[i], solstep.t̄[i], solstep.q̄[i], solstep.v̄[i])
end

function initialize!(solstep::SolutionStepPDAE, problem::AbstractProblemPDAE, extrap::Extrapolation = default_extrapolation())
    solstep.t  = initial_conditions(problem).t
    solstep.q .= initial_conditions(problem).q
    solstep.p .= initial_conditions(problem).p
    solstep.λ .= initial_conditions(problem).λ

    solstep.q̃ .= 0
    solstep.p̃ .= 0

    update_vector_fields!(solstep, problem)
    # update_implicit_functions!(solstep, problem)

    for i in eachhistory(solstep)
        solstep.t̄[i] = solstep.t̄[i-1] - timestep(problem)
        extrapolate!(solstep.t̄[i-1], solstep.q̄[i-1], solstep.p̄[i-1], solstep.t̄[i], solstep.q̄[i], solstep.p̄[i], problem, extrap)
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
        solstep.t̄[i]  = solstep.t̄[i-1]
        solstep.q̄[i] .= solstep.q̄[i-1]
        solstep.p̄[i] .= solstep.p̄[i-1]
        solstep.λ̄[i] .= solstep.λ̄[i-1]
        solstep.v̄[i] .= solstep.v̄[i-1]
        solstep.f̄[i] .= solstep.f̄[i-1]
        solstep.ū[i] .= solstep.ū[i-1]
        solstep.ḡ[i] .= solstep.ḡ[i-1]
    end

    solstep.t = solstep.t̄[0] += Δt

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
