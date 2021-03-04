@doc raw"""
`HODE`: Hamiltonian Ordinary Differential Equation *EXPERIMENTAL*

Defines a Hamiltonian ordinary differential initial value problem, that is
a canonical Hamiltonian system of equations,
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t), p(t)) , & q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), p(t)) , & p(t_{0}) &= p_{0} ,
\end{aligned}
```
with vector fields ``v`` and ``f``, given by
```math
\begin{aligned}
v &=   \frac{\partial H}{\partial p} , &
f &= - \frac{\partial H}{\partial q} ,
\end{aligned}
```
initial conditions ``(q_{0}, p_{0})`` and the dynamical variables ``(q,p)``
taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}``.

### Parameters

* `DT <: Number`: data type
* `TT <: Real`: time step type
* `AT <: AbstractArray{DT}`: array type
* `vType <: Function`: type of `v`
* `fType <: Function`: type of `f`
* `PType <: Function`: type of `P`
* `hamType <: Function`: Hamiltonian type
* `invType <: OptionalNamedTuple`: invariants type
* `parType <: OptionalNamedTuple`: parameters type
* `perType <: OptionalArray{AT}`: periodicity type

### Fields

* `d`: dimension of dynamical variables ``q`` and ``p`` as well as the vector fields ``v`` and ``f``
* `v`: function computing the vector field ``v``
* `f`: function computing the vector field ``f``
* `P`: function computing the Poisson matrix ``P``
* `t₀`: initial time
* `q₀`: initial condition for dynamical variable ``q``
* `p₀`: initial condition for dynamical variable ``p``
* `hamiltonian`: function computing the Hamiltonian ``H``
* `invariants`: either a `NamedTuple` containing the equation's invariants or `nothing`
* `parameters`: either a `NamedTuple` containing the equation's parameters or `nothing`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions

### Constructors

```julia
HODE(v, f, poisson, t₀, q₀, p₀, hamiltonian, invariants, parameters, periodicity)

HODE(v, f, h, t₀, q₀::StateVector, p₀::StateVector; kwargs...)
HODE(v, f, h, q₀::StateVector, p₀::StateVector; kwargs...)
HODE(v, f, h, t₀, q₀::State, p₀::State; kwargs...)
HODE(v, f, h, q₀::State, p₀::State; kwargs...)
```

### Keyword arguments:

* `poisson = symplectic_matrix`
* `invariants = nothing`
* `parameters = nothing`
* `periodicity = nothing`

"""
struct HODE{dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
            vType <: Function, fType <: Function, 
            PType <: Function,
            hamType <: Function,
            invType <: OptionalNamedTuple,
            parType <: OptionalNamedTuple,
            perType <: OptionalArray{arrayType}} <: AbstractEquationPODE{dType, tType}

    d::Int

    v::vType
    f::fType
    P::PType

    t₀::tType
    q₀::Vector{arrayType}
    p₀::Vector{arrayType}

    hamiltonian::hamType
    invariants::invType
    parameters::parType
    periodicity::perType

    function HODE(v, f, poisson, 
                  t₀::tType, q₀::Vector{arrayType}, p₀::Vector{arrayType},
                  hamiltonian, invariants, parameters, periodicity) where {
                        dType <: Number, tType <: Real, arrayType <: AbstractArray{dType}}

        d = length(q₀[begin])

        @assert length(q₀) == length(p₀)
        @assert all(length(q) == d for q in q₀)
        @assert all(length(p) == d for p in p₀)

        new{dType, tType, arrayType, typeof(v), typeof(f), typeof(poisson),
            typeof(hamiltonian), typeof(invariants), typeof(parameters), typeof(periodicity)}(
                d, v, f, poisson, t₀, q₀, p₀,
                hamiltonian, invariants, parameters, periodicity)
    end
end

_HODE(v, f, hamiltonian, t₀, q₀, p₀; poisson=symplectic_matrix, invariants=nothing, parameters=nothing, periodicity=nothing) = HODE(v, f, poisson, t₀, q₀, p₀, hamiltonian, invariants, parameters, periodicity)

HODE(v, f, h, t₀, q₀::StateVector, p₀::StateVector; kwargs...) = _HODE(v, f, h, t₀, q₀, p₀; kwargs...)
HODE(v, f, h, q₀::StateVector, p₀::StateVector; kwargs...) = HODE(v, f, h, 0.0, q₀, p₀; kwargs...)
HODE(v, f, h, t₀, q₀::State, p₀::State; kwargs...) = HODE(v, f, h, t₀, [q₀], [p₀]; kwargs...)
HODE(v, f, h, q₀::State, p₀::State; kwargs...) = HODE(v, f, h, 0.0, q₀, p₀; kwargs...)

const HODEinvType{invT,DT,TT,AT,VT,FT,PT,hamT,parT,perT} = HODE{DT,TT,AT,VT,FT,PT,hamT,invT,parT,perT} # type alias for dispatch on invariants type parameter
const HODEparType{parT,DT,TT,AT,VT,FT,PT,hamT,invT,perT} = HODE{DT,TT,AT,VT,FT,PT,hamT,invT,parT,perT} # type alias for dispatch on parameters type parameter
const HODEperType{perT,DT,TT,AT,VT,FT,PT,hamT,invT,parT} = HODE{DT,TT,AT,VT,FT,PT,hamT,invT,parT,perT} # type alias for dispatch on periodicity type parameter

Base.hash(ode::HODE, h::UInt) = hash(ode.d,
                        hash(ode.v, hash(ode.f, hash(ode.P,
                        hash(ode.t₀, hash(ode.q₀, hash(ode.p₀,
                        hash(ode.hamiltonian, hash(ode.invariants, hash(ode.parameters, hash(ode.periodicity, h)))))))))))

Base.:(==)(ode1::HODE, ode2::HODE) = (
                                ode1.d == ode2.d
                             && ode1.v == ode2.v
                             && ode1.f == ode2.f
                             && ode1.P == ode2.P
                             && ode1.t₀ == ode2.t₀
                             && ode1.q₀ == ode2.q₀
                             && ode1.p₀ == ode2.p₀
                             && ode1.hamiltonian == ode2.hamiltonian
                             && ode1.invariants  == ode2.invariants
                             && ode1.parameters  == ode2.parameters
                             && ode1.periodicity == ode2.periodicity)

function Base.similar(equ::HODE, t₀::Real, q₀::StateVector, p₀::StateVector;
                      parameters=equ.parameters)
    @assert all([length(q) == ndims(equ) for q in q₀])
    @assert all([length(p) == ndims(equ) for p in p₀])
    _HODE(equ.v, equ.f, equ.hamiltonian, t₀, q₀, p₀; poisson=equ.P, invariants=equ.invariants, parameters=parameters, periodicity=equ.periodicity)
end

Base.similar(equ::HODE, q₀, p₀; kwargs...) = similar(equ, equ.t₀, q₀, p₀; kwargs...)
Base.similar(equ::HODE, t₀::Real, q₀::State, p₀::State; kwargs...) = similar(equ, t₀, [q₀], [p₀]; kwargs...)

hashamiltonian(::HODE) = true

hasinvariants(::HODEinvType{<:Nothing}) = false
hasinvariants(::HODEinvType{<:NamedTuple}) = true

hasparameters(::HODEparType{<:Nothing}) = false
hasparameters(::HODEparType{<:NamedTuple}) = true

hasperiodicity(::HODEperType{<:Nothing}) = false
hasperiodicity(::HODEperType{<:AbstractArray}) = true

Base.axes(equ::HODE) = axes(equ.q₀[begin])
Base.ndims(equ::HODE) = equ.d
Common.nsamples(equ::HODE) = length(equ.q₀)

@inline Common.periodicity(equation::HODE) = hasperiodicity(equation) ? equation.periodicity : zero(equation.q₀[begin])
@inline initial_conditions(equation::HODE) = (equation.t₀, equation.q₀, equation.p₀)

_get_v(equ::HODE) = hasparameters(equ) ? (t,q,p,v) -> equ.v(t, q, p, v, equ.parameters) : equ.v
_get_f(equ::HODE) = hasparameters(equ) ? (t,q,p,f) -> equ.f(t, q, p, f, equ.parameters) : equ.f
_get_h(equ::HODE) = hasparameters(equ) ? (t,q,p)   -> equ.hamiltonian(t, q, p, equ.parameters) : equ.hamiltonian
_get_v̄(equ::HODE) = _get_v(equ)
_get_f̄(equ::HODE) = _get_f(equ)

function get_function_tuple(equ::HODE)
    names = (:v,:f,:h)
    equs  = (_get_v(equ), _get_f(equ), _get_h(equ))
    NamedTuple{names}(equs)
end
