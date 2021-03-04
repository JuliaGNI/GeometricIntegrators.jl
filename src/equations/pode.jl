@doc raw"""
`PODE`: Partitioned Ordinary Differential Equation

Defines a partitioned initial value problem
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t), p(t)) , &
q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), p(t)) , &
p(t_{0}) &= p_{0} ,
\end{aligned}
```
with vector fields ``v`` and ``f``, initial conditions ``(q_{0}, p_{0})`` and the solution
``(q,p)`` taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}``.

### Parameters

* `DT <: Number`: data type
* `TT <: Real`: time step type
* `AT <: AbstractArray{DT}`: array type
* `vType <: Function`: type of `v`
* `fType <: Function`: type of `f`
* `invType <: OptionalNamedTuple`: invariants type
* `parType <: OptionalNamedTuple`: parameters type
* `perType <: OptionalArray{AT}`: periodicity type

### Fields

* `d`: dimension of dynamical variables ``q`` and ``p`` as well as the vector fields ``v`` and ``f``
* `v`: function computing the vector field ``v``
* `f`: function computing the vector field ``f``
* `t₀`: initial time
* `q₀`: initial condition for `q`
* `p₀`: initial condition for `p`
* `invariants`: either a `NamedTuple` containing the equation's invariants or `nothing`
* `parameters`: either a `NamedTuple` containing the equation's parameters or `nothing`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions

The functions `v` and `f` must have the interface
```julia
    function v(t, q, p, v)
        v[1] = ...
        v[2] = ...
        ...
    end
```
and
```julia
    function f(t, q, p, f)
        f[1] = ...
        f[2] = ...
        ...
    end
```
where `t` is the current time, `q` and `p` are the current solution vectors
and `v` and `f` are the vectors which hold the result of evaluating the
vector fields ``v`` and ``f`` on `t`, `q` and `p`.

### Constructors

```julia
PODE(v, f, t₀, q₀, p₀, invariants, parameters, periodicity)

PODE(v, f, h, t₀, q₀::StateVector, p₀::StateVector; kwargs...)
PODE(v, f, q₀::StateVector, p₀::StateVector; kwargs...)
PODE(v, f, t₀, q₀::State, p₀::State; kwargs...)
PODE(v, f, q₀::State, p₀::State; kwargs...)
```

### Keyword arguments:

* `invariants = nothing`
* `parameters = nothing`
* `periodicity = nothing`

"""
struct PODE{dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
            vType <: Function, fType <: Function,
            invType <: OptionalNamedTuple,
            parType <: OptionalNamedTuple,
            perType <: OptionalArray{arrayType}} <: AbstractEquationPODE{dType, tType}

    d::Int

    v::vType
    f::fType

    t₀::tType
    q₀::Vector{arrayType}
    p₀::Vector{arrayType}

    invariants::invType
    parameters::parType
    periodicity::perType

    function PODE(v, f,
                  t₀::tType, q₀::Vector{arrayType}, p₀::Vector{arrayType},
                  invariants, parameters, periodicity) where {
                    dType <: Number, tType <: Real, arrayType <: AbstractArray{dType}}

        d = length(q₀[begin])

        @assert length(q₀) == length(p₀)
        @assert all(length(q) == d for q in q₀)
        @assert all(length(p) == d for p in p₀)

        new{dType, tType, arrayType,
            typeof(v), typeof(f),
            typeof(invariants), typeof(parameters), typeof(periodicity)}(
                d, v, f, t₀, q₀, p₀,
                invariants, parameters, periodicity)
    end
end

_PODE(v, f, t₀, q₀, p₀; invariants=nothing, parameters=nothing, periodicity=nothing) = PODE(v, f, t₀, q₀, p₀, invariants, parameters, periodicity)

PODE(v, f, t₀, q₀::StateVector, p₀::StateVector; kwargs...) = _PODE(v, f, t₀, q₀, p₀; kwargs...)
PODE(v, f, q₀::StateVector, p₀::StateVector; kwargs...) = PODE(v, f, 0.0, q₀, p₀; kwargs...)
PODE(v, f, t₀, q₀::State, p₀::State; kwargs...) = PODE(v, f, t₀, [q₀], [p₀]; kwargs...)
PODE(v, f, q₀::State, p₀::State; kwargs...) = PODE(v, f, 0.0, q₀, p₀; kwargs...)

const PODEinvType{invT,DT,TT,AT,VT,FT,parT,perT} = PODE{DT,TT,AT,VT,FT,invT,parT,perT} # type alias for dispatch on invariants type parameter
const PODEparType{parT,DT,TT,AT,VT,FT,invT,perT} = PODE{DT,TT,AT,VT,FT,invT,parT,perT} # type alias for dispatch on parameters type parameter
const PODEperType{perT,DT,TT,AT,VT,FT,invT,parT} = PODE{DT,TT,AT,VT,FT,invT,parT,perT} # type alias for dispatch on periodicity type parameter


Base.hash(ode::PODE, h::UInt) = hash(ode.d,
                        hash(ode.v, hash(ode.f,
                        hash(ode.t₀, hash(ode.q₀, hash(ode.p₀,
                        hash(ode.invariants, hash(ode.parameters, hash(ode.periodicity, h)))))))))

Base.:(==)(ode1::PODE, ode2::PODE) = (
                                ode1.d == ode2.d
                             && ode1.v == ode2.v
                             && ode1.f == ode2.f
                             && ode1.t₀ == ode2.t₀
                             && ode1.q₀ == ode2.q₀
                             && ode1.p₀ == ode2.p₀
                             && ode1.invariants  == ode2.invariants
                             && ode1.parameters  == ode2.parameters
                             && ode1.periodicity == ode2.periodicity)

function Base.similar(equ::PODE, t₀::Real, q₀::StateVector, p₀::StateVector; parameters=equ.parameters)
    @assert all([length(q) == ndims(equ) for q in q₀])
    @assert all([length(p) == ndims(equ) for p in p₀])
    PODE(equ.v, equ.f, t₀, q₀, p₀; invariants=equ.invariants, parameters=parameters, periodicity=equ.periodicity)
end

Base.similar(equ::PODE, q₀, p₀; kwargs...) = similar(equ, equ.t₀, q₀, p₀; kwargs...)
Base.similar(equ::PODE, t₀::Real, q₀::State, p₀::State; kwargs...) = similar(equ, t₀, [q₀], [p₀]; kwargs...)

hashamiltonian(::PODE) = false

hasinvariants(::PODEinvType{<:Nothing}) = false
hasinvariants(::PODEinvType{<:NamedTuple}) = true

hasparameters(::PODEparType{<:Nothing}) = false
hasparameters(::PODEparType{<:NamedTuple}) = true

hasperiodicity(::PODEperType{<:Nothing}) = false
hasperiodicity(::PODEperType{<:AbstractArray}) = true

Base.axes(equ::PODE) = axes(equ.q₀[begin])
Base.ndims(equ::PODE) = equ.d
Common.nsamples(equ::PODE) = length(equ.q₀)

@inline Common.periodicity(equation::PODE) = hasperiodicity(equation) ? equation.periodicity : zero(equation.q₀[begin])
@inline initial_conditions(equation::PODE) = (equation.t₀, equation.q₀, equation.p₀)

_get_v(equ::PODE) = hasparameters(equ) ? (t,q,p,v) -> equ.v(t, q, p, v, equ.parameters) : equ.v
_get_f(equ::PODE) = hasparameters(equ) ? (t,q,p,f) -> equ.f(t, q, p, f, equ.parameters) : equ.f
_get_v̄(equ::PODE) = _get_v(equ)
_get_f̄(equ::PODE) = _get_f(equ)


function get_functions(equ::PODE)
    names = (:v,:f)
    equs  = (_get_v(equ), _get_f(equ))

    NamedTuple{names}(equs)
end
