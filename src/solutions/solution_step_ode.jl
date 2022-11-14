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
* `q̃`: compensated summation error of q
* `v`: vector field of q
* `v̄`: vector field of q̄
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
            IT <: NamedTuple} <: SolutionStep{DT,TT,AT}

    t::TT
    t̄::TT

    q::AT
    q̄::AT
    q̃::AT

    v::VT
    v̄::VT

    internal::IT

    function SolutionStepODE(t::TT, q::AT, internal::IT = NamedTuple()) where {DT, TT, AT <: AbstractArray{DT}, IT}
        v = vectorfield(q)
        v̄ = vectorfield(q)
        new{DT,TT,AT,typeof(v),IT}(
            copy(t), zero(t),
            copy(q), zero(q), zero(q), 
            v, v̄, internal)
    end
end

function current(solstep::SolutionStepODE)
    (t = solstep.t, q = solstep.q)
end

function previous(solstep::SolutionStepODE)
    (t = solstep.t̄, q = solstep.q̄)
end

function Base.copy!(solstep::SolutionStepODE, sol::NamedTuple)
    solstep.t  = sol.t
    solstep.t̄  = 0
    solstep.q .= sol.q
    solstep.q̄ .= 0
    solstep.q̃ .= 0
    solstep.v .= 0
    solstep.v̄ .= 0
    return solstep
end

function GeometricBase.reset!(solstep::SolutionStepODE)
    solstep.t̄  = solstep.t
    solstep.q̄ .= solstep.q
    solstep.v̄ .= solstep.v
    return solstep
end

function update!(solstep::SolutionStepODE{DT,TT,AT}, Δt::TT, Δq::AT) where {DT,TT,AT}
    solstep.t += Δt
    for k in eachindex(Δq)
        solstep.q[k], solstep.q̃[k] = compensated_summation(Δq[k], solstep.q[k], solstep.q̃[k])
    end
    return solstep
end
