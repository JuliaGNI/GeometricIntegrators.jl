"""
Atomic solution for an DAE.

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
* `q̃`: compensated summation error of q
* `λ`: current solution of λ
* `λ̄`: previous solution of λ
* `v`: vector field of q
* `v̄`: vector field of q̄
* `u`: projective vector field of q
* `ū`: projective vector field of q̄
* `internal`: internal variables of the integrator (e.g., internal stages of a Runge-Kutta methods or solver output)

### Constructors

```julia
AtomicSolutionDAE{DT,TT,AT,IT}(nd, nm, internal::IT)
AtomicSolutionDAE{DT,TT,AT,IT}(t::TT, q::AT, λ::AT, internal::IT)
AtomicSolutionDAE(DT, TT, AT, nd, nm, internal::IT=NamedTuple())
AtomicSolutionDAE(t::TT, q::AT, λ::AT, internal::IT=NamedTuple())
```

* `nd`: dimension of the state vector
* `nm`: dimension of the constraint

"""
mutable struct AtomicSolutionDAE{DT <: Number, TT <: Real, AT <: AbstractArray{DT}, IT <: NamedTuple} <: AtomicSolution{DT,TT,AT}
    t::TT
    t̄::TT

    q::AT
    q̄::AT
    q̃::AT

    λ::AT
    λ̄::AT

    v::AT
    v̄::AT

    u::AT
    ū::AT

    internal::IT

    function AtomicSolutionDAE{DT,TT,AT,IT}(nd, nm, internal::IT) where {DT,TT,AT,IT}
        new(zero(TT), zero(TT),
            AT(zeros(DT, nd)), AT(zeros(DT, nd)), AT(zeros(DT, nd)),
            AT(zeros(DT, nm)), AT(zeros(DT, nm)),
            AT(zeros(DT, nd)), AT(zeros(DT, nd)),
            AT(zeros(DT, nd)), AT(zeros(DT, nd)),
            internal)
    end

    function AtomicSolutionDAE{DT,TT,AT,IT}(t::TT, q::AT, λ::AT, internal::IT) where {DT,TT,AT,IT}
        new(zero(t), zero(t),
            zero(q), zero(q), zero(q),
            zero(λ), zero(λ),
            zero(q), zero(q),
            zero(λ), zero(λ),
            internal)
    end
end

AtomicSolutionDAE(DT, TT, AT, nd, nm, internal::IT=NamedTuple()) where {IT} = AtomicSolutionDAE{DT,TT,AT,IT}(nd, nm, internal)
AtomicSolutionDAE(t::TT, q::AT, λ::AT, internal::IT=NamedTuple()) where {DT, TT, AT <: AbstractArray{DT}, IT} = AtomicSolutionDAE{DT,TT,AT,IT}(t, q, λ, internal)

function set_initial_conditions!(asol::AtomicSolutionDAE, equ::AbstractEquationDAE, i::Int=1)
    @assert i ≥ nsamples(equ)
    t, q, λ = initial_conditions(equ)
    asol.t  = t
    asol.q .= q[i]
    asol.λ .= λ[i]
    asol.v .= 0
end

function set_solution!(asol::AtomicSolutionDAE, sol)
    t, q, λ = sol
    asol.t  = t
    asol.q .= q
    asol.λ .= λ
    asol.v .= 0
end

function get_solution(asol::AtomicSolutionDAE)
    (asol.t, asol.q, asol.λ)
end

function GeometricBase.reset!(asol::AtomicSolutionDAE, Δt)
    asol.t̄  = asol.t
    asol.q̄ .= asol.q
    asol.λ̄ .= asol.λ
    asol.v̄ .= asol.v
    asol.ū .= asol.u
    asol.t += Δt
end

function update!(asol::AtomicSolutionDAE{DT}, y::Vector{DT}, λ::Vector{DT}) where {DT}
    for k in eachindex(y)
        update!(asol, y[k], k)
    end
    for k in eachindex(λ)
        asol.λ[k] = λ[k]
    end
end

function update!(asol::AtomicSolutionDAE{DT}, y::DT, k::Union{Int,CartesianIndex}) where {DT}
    asol.q[k], asol.q̃[k] = compensated_summation(y, asol.q[k], asol.q̃[k])
end
