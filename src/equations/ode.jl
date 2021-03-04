@doc raw"""
`ODE`: Ordinary Differential Equation

Defines an initial value problem
```math
\dot{q} (t) = v(t, q(t)) , \qquad q(t_{0}) = q_{0} ,
```
with vector field ``v``, initial condition ``q_{0}`` and the solution
``q`` taking values in ``\mathbb{R}^{d}``.

### Parameters

* `DT <: Number`: data type
* `TT <: Real`: time step type
* `AT <: AbstractArray{DT}`: array type
* `vType <: Function`: type of `v`
* `invType <: OptionalNamedTuple`: invariants type
* `parType <: OptionalNamedTuple`: parameters type
* `perType <: OptionalArray{AT}`: periodicity type

### Fields

* `d`: dimension of dynamical variable ``q`` and the vector field ``v``
* `v`: function computing the vector field
* `t₀`: initial time
* `q₀`: initial condition
* `invariants`: either a `NamedTuple` containing the equation's invariants or `nothing`
* `parameters`: either a `NamedTuple` containing the equation's parameters or `nothing`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions

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

### Constructors

```julia
ODE(v, t₀, q₀, invariants, parameters, periodicity)

ODE(v, t₀, q₀::StateVector; kwargs...)
ODE(v, q₀::StateVector; kwargs...)
ODE(v, t₀, q₀::State; kwargs...)
ODE(v, q₀::State; kwargs...)
```

### Keyword arguments:

* `invariants = nothing`
* `parameters = nothing`
* `periodicity = nothing`

"""
struct ODE{dType <: Number, tType <: Real, arrayType <: AbstractArray{dType}, vType <: Function,
            invType <: OptionalNamedTuple,
            parType <: OptionalNamedTuple,
            perType <: OptionalArray{arrayType}} <: AbstractEquationODE{dType, tType}

    d::Int

    v::vType

    t₀::tType
    q₀::Vector{arrayType}

    invariants::invType
    parameters::parType
    periodicity::perType

    function ODE(v, t₀::tType, q₀::Vector{arrayType},
                 invariants, parameters, periodicity) where {
                        dType <: Number, tType <: Real, arrayType <: AbstractArray{dType}}

        d = length(q₀[begin])

        @assert all(length(q) == d for q in q₀)

        new{dType, tType, arrayType, typeof(v),
            typeof(invariants), typeof(parameters), typeof(periodicity)}(
                d, v, t₀, q₀, invariants, parameters, periodicity)
    end
end

_ODE(v, t₀, q₀; invariants=nothing, parameters=nothing, periodicity=nothing) = ODE(v, t₀, q₀, invariants, parameters, periodicity)

ODE(v, t₀, q₀::StateVector; kwargs...) = _ODE(v, t₀, q₀; kwargs...)
ODE(v, q₀::StateVector; kwargs...) = ODE(v, 0.0, q₀; kwargs...)
ODE(v, t₀, q₀::State; kwargs...) = ODE(v, t₀, [q₀]; kwargs...)
ODE(v, q₀::State; kwargs...) = ODE(v, 0.0, q₀; kwargs...)

const ODEinvType{invT,DT,TT,AT,VT,parT,perT} = ODE{DT,TT,AT,VT,invT,parT,perT} # type alias for dispatch on invariants type parameter
const ODEparType{parT,DT,TT,AT,VT,invT,perT} = ODE{DT,TT,AT,VT,invT,parT,perT} # type alias for dispatch on parameters type parameter
const ODEperType{perT,DT,TT,AT,VT,invT,parT} = ODE{DT,TT,AT,VT,invT,parT,perT} # type alias for dispatch on periodicity type parameter


Base.hash(ode::ODE, h::UInt) = hash(ode.d, hash(ode.v, hash(ode.t₀, hash(ode.q₀, 
                        hash(ode.invariants, hash(ode.parameters, hash(ode.periodicity, h)))))))

Base.:(==)(ode1::ODE, ode2::ODE) = (
                                ode1.d == ode2.d
                             && ode1.v == ode2.v
                             && ode1.t₀ == ode2.t₀
                             && ode1.q₀ == ode2.q₀
                             && ode1.invariants  == ode2.invariants
                             && ode1.parameters  == ode2.parameters
                             && ode1.periodicity == ode2.periodicity)

function Base.similar(equ::ODE, t₀::Real, q₀::StateVector; parameters=equ.parameters)
    @assert all([length(q) == ndims(equ) for q in q₀])
    ODE(equ.v, t₀, q₀; invariants=equ.invariants, parameters=parameters, periodicity=equ.periodicity)
end

Base.similar(equ::ODE, q₀; kwargs...) = similar(equ, equ.t₀, q₀; kwargs...)
Base.similar(equ::ODE, t₀::Real, q₀::State; kwargs...) = similar(equ, t₀, [q₀]; kwargs...)

hashamiltonian(::ODE) = false

hasinvariants(::ODEinvType{<:Nothing}) = false
hasinvariants(::ODEinvType{<:NamedTuple}) = true

hasparameters(::ODEparType{<:Nothing}) = false
hasparameters(::ODEparType{<:NamedTuple}) = true

hasperiodicity(::ODEperType{<:Nothing}) = false
hasperiodicity(::ODEperType{<:AbstractArray}) = true

Base.axes(equ::ODE) = axes(equ.q₀[begin])
Base.ndims(equ::ODE) = equ.d
Common.nsamples(equ::ODE) = length(equ.q₀)

@inline Common.periodicity(equation::ODE) = hasperiodicity(equation) ? equation.periodicity : zero(equation.q₀[begin])
@inline initial_conditions(equation::ODE) = (equation.t₀, equation.q₀)

_get_v(equ::ODE) = hasparameters(equ) ? (t,q,v) -> equ.v(t, q, v, equ.parameters) : equ.v
_get_v̄(equ::ODE) = _get_v(equ)

function get_function_tuple(equ::ODE)
    names = (:v,)
    equs  = (_get_v(equ),)

    NamedTuple{names}(equs)
end

function get_invariants(equ::ODE)
    if hasinvariants(equ)
        keys = ()
        invs = ()
        for (key, inv) in pairs(equ.invariants)
            keys = (keys..., key)
            invs = (invs..., hasparameters(equ) ? (t,q) -> inv(t, q, equ.parameters) : inv)
        end
        return NamedTuple{keys}(invs)
    else
        return NamedTuple()
    end
end