"""
Atomic solution for an PODE.

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
* `v`: vector field of q
* `v̄`: vector field of q̄
* `f`: vector field of p
* `f̄`: vector field of p̄
* `internal`: internal variables of the integrator (e.g., internal stages of a Runge-Kutta methods or solver output)

### Constructors

```julia
AtomicSolutionPODE{DT,TT,AT,IT}(nd, internal::IT)
AtomicSolutionPODE{DT,TT,AT,IT}(t::TT, q::AT, p::AT, internal::IT)
AtomicSolutionPODE(DT, TT, AT, nd, internal::IT=NamedTuple())
AtomicSolutionPODE(t::TT, q::AT, p::AT, internal::IT=NamedTuple())
```

* `nd`: dimension of the state vector

"""
mutable struct AtomicSolutionPODE{DT <: Number, TT <: Real, AT <: AbstractArray{DT}, IT <: NamedTuple} <: AtomicSolution{DT,TT,AT}
    t::TT
    t̄::TT

    q::AT
    q̄::AT
    q̃::AT

    p::AT
    p̄::AT
    p̃::AT

    v::AT
    v̄::AT
    f::AT
    f̄::AT

    internal::IT

    function AtomicSolutionPODE{DT,TT,AT,IT}(nd, internal::IT) where {DT,TT,AT,IT}
        new(zero(TT), zero(TT),
            AT(zeros(DT, nd)), AT(zeros(DT, nd)), AT(zeros(DT, nd)),
            AT(zeros(DT, nd)), AT(zeros(DT, nd)), AT(zeros(DT, nd)),
            AT(zeros(DT, nd)), AT(zeros(DT, nd)), AT(zeros(DT, nd)), AT(zeros(DT, nd)),
            internal)
    end

    function AtomicSolutionPODE{DT,TT,AT,IT}(t::TT, q::AT, p::AT, internal::IT) where {DT,TT,AT,IT}
        new(zero(t), zero(t),
            zero(q), zero(q), zero(q),
            zero(p), zero(p), zero(p),
            zero(q), zero(q), zero(p), zero(p),
            internal)
    end
end

AtomicSolutionPODE(DT, TT, AT, nd, internal::IT=NamedTuple()) where {IT} = AtomicSolutionPODE{DT,TT,AT,IT}(nd, internal)
AtomicSolutionPODE(t::TT, q::AT, p::AT, internal::IT=NamedTuple()) where {DT, TT, AT <: AbstractArray{DT}, IT} = AtomicSolutionPODE{DT,TT,AT,IT}(t, q, p, internal)

function set_initial_conditions!(asol::AtomicSolutionPODE, equ::AbstractEquationPODE, i::Int=1)
    @assert i ≥ nsamples(equ)
    t, q, p = initial_conditions(equ)
    asol.t  = t
    asol.q .= q[i]
    asol.p .= p[i]
    asol.v .= 0
    asol.f .= 0
end

function set_solution!(asol::AtomicSolutionPODE, sol)
    t, q, p = sol
    asol.t  = t
    asol.q .= q
    asol.p .= p
    asol.v .= 0
    asol.f .= 0
end

function get_solution(asol::AtomicSolutionPODE)
    (asol.t, asol.q, asol.p)
end

function GeometricBase.reset!(asol::AtomicSolutionPODE, Δt)
    asol.t̄  = asol.t
    asol.q̄ .= asol.q
    asol.p̄ .= asol.p
    asol.v̄ .= asol.v
    asol.f̄ .= asol.f
    asol.t += Δt
end

function update!(asol::AtomicSolutionPODE{DT}, y::Vector{DT}, z::Vector{DT}) where {DT}
    for k in eachindex(y,z)
        update!(asol, y[k], z[k], k)
    end
end

function update!(asol::AtomicSolutionPODE{DT}, y::DT, z::DT, k::Union{Int,CartesianIndex}) where {DT}
    asol.q[k], asol.q̃[k] = compensated_summation(y, asol.q[k], asol.q̃[k])
    asol.p[k], asol.p̃[k] = compensated_summation(z, asol.p[k], asol.p̃[k])
end
