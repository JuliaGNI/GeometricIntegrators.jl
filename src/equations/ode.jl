@doc raw"""
`ODE`: Ordinary Differential Equation

Defines an initial value problem
```math
\dot{q} (t) = v(t, q(t)) , \qquad q(t_{0}) = q_{0} ,
```
with vector field ``v``, initial condition ``q_{0}`` and the solution
``q`` taking values in ``\mathbb{R}^{d}``.

### Fields

* `d`: dimension of dynamical variable ``q`` and the vector field ``v``
* `v`: function computing the vector field
* `h`: function computing the Hamiltonian (optional)
* `t₀`: initial time
* `q₀`: initial condition

The function `v` providing the vector field must have the interface
```julia
    function v(t, q, v)
        v[1] = ...
        v[2] = ...
        ...
    end
```
where `t` is the current time, `q` is the current solution vector, and
`v` is the vector which holds the result of evaluating the vector field ``v``
on `t` and `q`.

"""
struct ODE{dType <: Number, tType <: Real, arrayType <: AbstractArray{dType}, vType <: Function,
           hType <: OptionalFunction, pType <: Union{NamedTuple,Nothing}} <: AbstractEquationODE{dType, tType}

    d::Int
    v::vType
    h::hType
    t₀::tType
    q₀::Vector{arrayType}
    parameters::pType
    periodicity::Vector{dType}

    function ODE(v::vType, t₀::tType, q₀::Vector{arrayType};
                 h::hType=nothing, parameters::pType=nothing,
                 periodicity=zero(q₀[begin])) where {
                        dType <: Number, tType <: Real, arrayType <: AbstractArray{dType}, vType <: Function,
                        hType <: OptionalFunction, pType <: Union{NamedTuple,Nothing}}

        d = length(q₀[begin])

        @assert all([length(q) == d for q in q₀])

        new{dType, tType, arrayType, vType, hType, pType}(d, v, h, t₀, q₀, parameters, periodicity)
    end
end

ODE(v, q₀::StateVector; kwargs...) = ODE(v, 0.0, q₀; kwargs...)
ODE(v, t₀, q₀::State; kwargs...) = ODE(v, t₀, [q₀]; kwargs...)
ODE(v, q₀::State; kwargs...) = ODE(v, 0.0, q₀; kwargs...)

const ODEHT{HT,DT,TT,AT,VT,PT} = ODE{DT,TT,AT,VT,HT,PT} # type alias for dispatch on Hamiltonian type parameter
const ODEPT{PT,DT,TT,AT,VT,HT} = ODE{DT,TT,AT,VT,HT,PT} # type alias for dispatch on parameters type parameter


Base.hash(ode::ODE, h::UInt) = hash(ode.d, hash(ode.v, hash(ode.t₀,
        hash(ode.q₀, hash(ode.periodicity, hash(ode.parameters, h))))))

Base.:(==)(ode1::ODE, ode2::ODE) = (
                                ode1.d == ode2.d
                             && ode1.v == ode2.v
                             && ode1.h == ode2.h
                             && ode1.t₀ == ode2.t₀
                             && ode1.q₀ == ode2.q₀
                             && ode1.parameters == ode2.parameters
                             && ode1.periodicity == ode2.periodicity)

Base.similar(equ::ODE, q₀; kwargs...) = similar(equ, equ.t₀, q₀; kwargs...)
Base.similar(equ::ODE, t₀::Real, q₀::State; kwargs...) = similar(equ, t₀, [q₀]; kwargs...)

function Base.similar(equ::ODE, t₀::Real, q₀::StateVector;
                      h=equ.h, parameters=equ.parameters, periodicity=equ.periodicity)
    @assert all([length(q) == ndims(equ) for q in q₀])
    ODE(equ.v, t₀, q₀; h=h, parameters=parameters, periodicity=periodicity)
end

Base.ndims(equ::ODE) = equ.d
Common.nsamples(equ::ODE) = length(equ.q₀)
Common.periodicity(equ::ODE) = equ.periodicity

initial_conditions(equ::ODE) = (equ.t₀, equ.q₀)

hashamiltonian(::ODEHT{<:Nothing}) = false
hashamiltonian(::ODEHT{<:Function}) = true

hasparameters(::ODEPT{<:Nothing}) = false
hasparameters(::ODEPT{<:NamedTuple}) = true

_get_v(equ::ODE) = hasparameters(equ) ? (t,q,v) -> equ.v(t, q, v, equ.parameters) : equ.v
_get_h(equ::ODE) = hasparameters(equ) ? (t,q) -> equ.h(t, q, equ.parameters) : equ.h


function get_function_tuple(equ::ODE)
    names = (:v,)
    equs  = (_get_v(equ),)

    if hashamiltonian(equ)
        names = (names..., :h)
        equs  = (equs..., _get_h(equ))
    end

    NamedTuple{names}(equs)
end
