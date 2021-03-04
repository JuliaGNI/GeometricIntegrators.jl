"""
Atomic solution for an PDAE.

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
* `p`: current solution of p
* `p̄`: previous solution of p
* `p̃`: compensated summation error of p
* `λ`: current solution of λ
* `λ̄`: previous solution of λ
* `v`: vector field of q
* `v̄`: vector field of q̄
* `f`: vector field of p
* `f̄`: vector field of p̄
* `u`: projective vector field of q
* `ū`: projective vector field of q̄
* `g`: projective vector field of p
* `ḡ`: projective vector field of p̄
* `internal`: internal variables of the integrator (e.g., internal stages of a Runge-Kutta methods or solver output)

### Constructors

```julia
AtomicSolutionPDAE{DT,TT,AT,IT}(nd, nm, internal::IT)
AtomicSolutionPDAE{DT,TT,AT,IT}(t::TT, q::AT, p::AT, λ::AT, internal::IT)
AtomicSolutionPDAE(DT, TT, AT, nd, nm, internal::IT=NamedTuple())
AtomicSolutionPDAE(t::TT, q::AT, p::AT, λ::AT, internal::IT=NamedTuple())
```

* `nd`: dimension of the state vector
* `nm`: dimension of the constraint

"""
mutable struct AtomicSolutionPDAE{DT <: Number, TT <: Real, AT <: AbstractArray{DT}, IT <: NamedTuple} <: AtomicSolution{DT,TT,AT}
    t::TT
    t̄::TT

    q::AT
    q̄::AT
    q̃::AT

    p::AT
    p̄::AT
    p̃::AT

    λ::AT
    λ̄::AT

    v::AT
    v̄::AT
    f::AT
    f̄::AT
    
    u::AT
    ū::AT
    g::AT
    ḡ::AT

    internal::IT

    function AtomicSolutionPDAE{DT,TT,AT,IT}(nd, nm, internal::IT) where {DT,TT,AT,IT}
        new(zero(TT), zero(TT),
            AT(zeros(DT, nd)), AT(zeros(DT, nd)), AT(zeros(DT, nd)),
            AT(zeros(DT, nd)), AT(zeros(DT, nd)), AT(zeros(DT, nd)),
            AT(zeros(DT, nm)), AT(zeros(DT, nm)),
            AT(zeros(DT, nd)), AT(zeros(DT, nd)), AT(zeros(DT, nd)), AT(zeros(DT, nd)),
            AT(zeros(DT, nd)), AT(zeros(DT, nd)), AT(zeros(DT, nd)), AT(zeros(DT, nd)),
            internal)
    end

    function AtomicSolutionPDAE{DT,TT,AT,IT}(t::TT, q::AT, p::AT, λ::AT, internal::IT) where {DT,TT,AT,IT}
        new(zero(t), zero(t),
            zero(q), zero(q), zero(q),
            zero(p), zero(p), zero(p),
            zero(λ), zero(λ),
            zero(q), zero(q), zero(p), zero(p),
            zero(λ), zero(λ), zero(λ), zero(λ),
            internal)
    end
end

AtomicSolutionPDAE(DT, TT, AT, nd, nm, internal::IT=NamedTuple()) where {IT} = AtomicSolutionPDAE{DT,TT,AT,IT}(nd, nm, internal)
AtomicSolutionPDAE(t::TT, q::AT, p::AT, λ::AT, internal::IT=NamedTuple()) where {DT, TT, AT <: AbstractArray{DT}, IT} = AtomicSolutionPDAE{DT,TT,AT,IT}(t, q, p, λ, internal)

function set_initial_conditions!(asol::AtomicSolutionPDAE, equ::AbstractEquationPDAE, i::Int=1)
    @assert i ≥ nsamples(equ)
    t, q, p, λ = initial_conditions(equ)
    asol.t  = t
    asol.q .= q[i]
    asol.p .= p[i]
    asol.λ .= λ[i]
    asol.v .= 0
    asol.f .= 0
end

function set_solution!(asol::AtomicSolutionPDAE, sol)
    t, q, p, λ = sol
    asol.t  = t
    asol.q .= q
    asol.p .= p
    asol.λ .= λ
    asol.v .= 0
    asol.f .= 0
end

function get_solution(asol::AtomicSolutionPDAE)
    (asol.t, asol.q, asol.p, asol.λ)
end

function Common.reset!(asol::AtomicSolutionPDAE, Δt)
    asol.t̄  = asol.t
    asol.q̄ .= asol.q
    asol.p̄ .= asol.p
    asol.λ̄.= asol.λ
    asol.v̄ .= asol.v
    asol.f̄ .= asol.f
    asol.ū .= asol.u
    asol.ḡ .= asol.g
    asol.t += Δt
end

function update!(asol::AtomicSolutionPDAE{DT}, y::Vector{DT}, z::Vector{DT}, λ::Vector{DT}) where {DT}
    for k in eachindex(y,z)
        update!(asol, y[k], z[k], k)
    end
    for k in eachindex(λ)
        asol.λ[k] = λ[k]
    end
end

function update!(asol::AtomicSolutionPDAE{DT}, y::DT, z::DT, k::Union{Int,CartesianIndex}) where {DT}
    asol.q[k], asol.q̃[k] = compensated_summation(y, asol.q[k], asol.q̃[k])
    asol.p[k], asol.p̃[k] = compensated_summation(z, asol.p[k], asol.p̃[k])
end
