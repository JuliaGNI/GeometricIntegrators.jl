@doc raw"""
`LODE`: Variational Ordinary Differential Equation *EXPERIMENTAL*

Defines an implicit initial value problem
```math
\begin{aligned}
\dot{q} (t) &= v(t) , &
q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), v(t)) , &
p(t_{0}) &= p_{0} , \\
p(t) &= ϑ(t, q(t), v(t))
\end{aligned}
```
with vector field ``f``, the momentum defined by ``p``, initial conditions ``(q_{0}, p_{0})`` and the solution
``(q,p)`` taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}``.
This is a special case of a differential algebraic equation with dynamical
variables ``(q,p)`` and algebraic variable ``v``.

### Parameters

* `DT <: Number`: data type
* `TT <: Real`: time step type
* `AT <: AbstractArray{DT}`: array type
* `ϑType <: Function`: type of `ϑ`
* `fType <: Function`: type of `f`
* `gType <: Function`: type of `g`
* `ωType <: Function`: type of `ω`
* `v̄Type <: Function`: type of `v̄`
* `f̄Type <: Function`: type of `f̄`
* `lagType <: Function`: Lagrangian type
* `invType <: OptionalNamedTuple`: invariants type
* `parType <: OptionalNamedTuple`: parameters type
* `perType <: OptionalArray{AT}`: periodicity type

### Fields

* `d`: dimension of dynamical variables ``q`` and ``p`` as well as the vector fields ``f`` and ``p``
* `ϑ`: function determining the momentum
* `f`: function computing the vector field
* `g`: function determining the projection, given by ``\nabla \vartheta (q) \cdot \lambda``
* `ω`: function computing the symplectic matrix
* `v̄`: function computing an initial guess for the velocity field ``v`` (optional)
* `f̄`: function computing an initial guess for the force field ``f`` (optional)
* `t₀`: initial time (optional)
* `q₀`: initial condition for `q`
* `p₀`: initial condition for `p`
* `λ₀`: initial condition for `λ` (optional)
* `lagrangian`: function computing the Lagrangian ``L``
* `invariants`: either a `NamedTuple` containing the equation's invariants or `nothing`
* `parameters`: either a `NamedTuple` containing the equation's parameters or `nothing`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions

The functions `ϑ` and `f` must have the interface
```julia
    function ϑ(t, q, v, p)
        p[1] = ...
        p[2] = ...
        ...
    end
```
and
```julia
    function f(t, q, v, f)
        f[1] = ...
        f[2] = ...
        ...
    end
```
where `t` is the current time, `q` is the current solution vector, `v` is the
current velocity and `f` and `p` are the vectors which hold the result of
evaluating the functions ``f`` and ``ϑ`` on `t`, `q` and `v`.
The funtions `g` and `v` are specified by
```julia
    function g(t, q, λ, g)
        g[1] = ...
        g[2] = ...
        ...
    end
```
and
```julia
    function v(t, q, p, v)
        v[1] = ...
        v[2] = ...
        ...
    end
```

### Constructors

```julia
LODE(ϑ, f, ω, v̄, f̄, t₀, q₀, p₀, λ₀, lagrangian, invariants, parameters, periodicity)

LODE(ϑ, f, l, ω, t₀, q₀::StateVector, p₀::StateVector, λ₀::StateVector=zero(q₀); kwargs...)
LODE(ϑ, f, l, ω, q₀::StateVector, p₀::StateVector, λ₀::StateVector=zero(q₀); kwargs...)
LODE(ϑ, f, l, ω, t₀, q₀::State, p₀::State, λ₀::StateVector=zero(q₀); kwargs...)
LODE(ϑ, f, l, ω, q₀::State, p₀::State, λ₀::StateVector=zero(q₀); kwargs...)
```

### Keyword arguments

* `v̄ = (t,q,v) -> nothing`
* `f̄ = f`
* `invariants = nothing`
* `parameters = nothing`
* `periodicity = nothing`

"""
struct LODE{dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
            ϑType <: Function, fType <: Function, gType <: Function,
            ωType <: Function, v̄Type <: Function, f̄Type <: Function,
            lagType <: Function,
            invType <: OptionalNamedTuple,
            parType <: OptionalNamedTuple,
            perType <: OptionalArray{arrayType}} <: AbstractEquationPODE{dType, tType}

    d::Int
    m::Int

    ϑ::ϑType
    f::fType
    g::gType
    ω::ωType
    v̄::v̄Type
    f̄::f̄Type

    t₀::tType
    q₀::Vector{arrayType}
    p₀::Vector{arrayType}
    λ₀::Vector{arrayType}

    lagrangian::lagType
    invariants::invType
    parameters::parType
    periodicity::perType

    function LODE(ϑ, f, g, ω, v̄, f̄,
                t₀::tType, q₀::Vector{arrayType}, p₀::Vector{arrayType}, λ₀::Vector{arrayType},
                lagrangian, invariants, parameters, periodicity) where {
                    dType <: Number, tType <: Real, arrayType <: AbstractArray{dType}}

        d = length(q₀[begin])

        @assert length(q₀) == length(p₀)
        @assert all(length(q) == d for q in q₀)
        @assert all(length(p) == d for p in p₀)
        @assert all(length(λ) == d for λ in λ₀)

        new{dType, tType, arrayType,
            typeof(ϑ), typeof(f), typeof(g), typeof(ω), typeof(v̄), typeof(f̄),
            typeof(lagrangian), typeof(invariants), typeof(parameters), typeof(periodicity)}(
                d, d, ϑ, f, g, ω, v̄, f̄, t₀, q₀, p₀, λ₀, lagrangian, invariants, parameters, periodicity)
    end
end

_LODE(ϑ, f, g, lagrangian, ω, t₀, q₀, p₀, λ₀; v̄=(parameters === nothing ? (t,q,v)->nothing :  (t,q,v,params)->nothing), f̄=f, invariants=nothing, parameters=nothing, periodicity=nothing) = LODE(ϑ, f, g, ω, v̄, f̄, t₀, q₀, p₀, λ₀, lagrangian, invariants, parameters, periodicity)

LODE(ϑ, f, g, l, ω, t₀, q₀::StateVector, p₀::StateVector, λ₀::StateVector=zero(q₀); kwargs...) = _LODE(ϑ, f, g, l, ω, t₀, q₀, p₀, λ₀; kwargs...)
LODE(ϑ, f, g, l, ω, q₀::StateVector, p₀::StateVector, λ₀::StateVector=zero(q₀); kwargs...) = LODE(ϑ, f, g, l, ω, 0.0, q₀, p₀, λ₀; kwargs...)
LODE(ϑ, f, g, l, ω, t₀, q₀::State, p₀::State, λ₀::State=zero(q₀); kwargs...) = LODE(ϑ, f, g, l, ω, t₀, [q₀], [p₀], [λ₀]; kwargs...)
LODE(ϑ, f, g, l, ω, q₀::State, p₀::State, λ₀::State=zero(q₀); kwargs...) = LODE(ϑ, f, g, l, ω, 0.0, q₀, p₀, λ₀; kwargs...)

const LODEinvType{invT,DT,TT,AT,ΘT,FT,GT,ΩT,ŪT,ḠT,lagT,parT,perT} = LODE{DT,TT,AT,ΘT,FT,GT,ΩT,ŪT,ḠT,lagT,invT,parT,perT} # type alias for dispatch on invariants type parameter
const LODEparType{parT,DT,TT,AT,ΘT,FT,GT,ΩT,ŪT,ḠT,lagT,invT,perT} = LODE{DT,TT,AT,ΘT,FT,GT,ΩT,ŪT,ḠT,lagT,invT,parT,perT} # type alias for dispatch on parameters type parameter
const LODEperType{perT,DT,TT,AT,ΘT,FT,GT,ΩT,ŪT,ḠT,lagT,invT,parT} = LODE{DT,TT,AT,ΘT,FT,GT,ΩT,ŪT,ḠT,lagT,invT,parT,perT} # type alias for dispatch on periodicity type parameter

Base.hash(ode::LODE, h::UInt) = hash(ode.d, hash(ode.m,
          hash(ode.ϑ, hash(ode.f, hash(ode.g,
          hash(ode.ω, hash(ode.v̄, hash(ode.f̄,
          hash(ode.t₀, hash(ode.q₀, hash(ode.p₀, hash(ode.λ₀,
          hash(ode.invariants, hash(ode.parameters, hash(ode.periodicity, h)))))))))))))))

Base.:(==)(ode1::LODE, ode2::LODE) = (
                                ode1.d == ode2.d
                             && ode1.m == ode2.m
                             && ode1.ϑ == ode2.ϑ
                             && ode1.f == ode2.f
                             && ode1.g == ode2.g
                             && ode1.ω == ode2.ω
                             && ode1.v̄ == ode2.v̄
                             && ode1.f̄ == ode2.f̄
                             && ode1.t₀ == ode2.t₀
                             && ode1.q₀ == ode2.q₀
                             && ode1.p₀ == ode2.p₀
                             && ode1.λ₀ == ode2.λ₀
                             && ode1.lagrangian == ode2.lagrangian
                             && ode1.invariants  == ode2.invariants
                             && ode1.parameters  == ode2.parameters
                             && ode1.periodicity == ode2.periodicity)

function Base.similar(equ::LODE, t₀::Real, q₀::StateVector, p₀::StateVector, λ₀::StateVector; parameters=equ.parameters)
    @assert all([length(q) == ndims(equ) for q in q₀])
    @assert all([length(p) == ndims(equ) for p in p₀])
    @assert all([length(λ) == ndims(equ) for λ in λ₀])
    LODE(equ.ϑ, equ.f, equ.g, equ.lagrangian, equ.ω, t₀, q₀, p₀, λ₀; v̄=equ.v̄, f̄=equ.f̄, invariants=equ.invariants, parameters=parameters, periodicity=equ.periodicity)
end

Base.similar(equ::LODE, q₀, p₀, λ₀=get_λ₀(q₀, equ.λ₀); kwargs...) = similar(equ, equ.t₀, q₀, p₀, λ₀; kwargs...)
Base.similar(equ::LODE, t₀::Real, q₀::State, p₀::State, λ₀::State=get_λ₀(q₀, equ.λ₀); kwargs...) = similar(equ, t₀, [q₀], [p₀], [λ₀]; kwargs...)

hasinvariants(::LODEinvType{<:Nothing}) = false
hasinvariants(::LODEinvType{<:NamedTuple}) = true

hasparameters(::LODEparType{<:Nothing}) = false
hasparameters(::LODEparType{<:NamedTuple}) = true

hasperiodicity(::LODEperType{<:Nothing}) = false
hasperiodicity(::LODEperType{<:AbstractArray}) = true

Base.axes(equ::LODE) = axes(equ.q₀[begin])
Base.ndims(equ::LODE) = equ.d
Common.nsamples(equ::LODE) = length(equ.q₀)

@inline Common.periodicity(equation::LODE) = hasperiodicity(equation) ? equation.periodicity : zero(equation.q₀[begin])
@inline initial_conditions(equation::LODE) = (equation.t₀, equation.q₀, equation.p₀, equation.λ₀)

_get_ϑ(equ::LODE) = hasparameters(equ) ? (t,q,v,ϑ) -> equ.ϑ(t, q, v, ϑ, equ.parameters) : equ.ϑ
_get_f(equ::LODE) = hasparameters(equ) ? (t,q,v,f) -> equ.f(t, q, v, f, equ.parameters) : equ.f
_get_g(equ::LODE) = hasparameters(equ) ? (t,q,v,g) -> equ.g(t, q, v, g, equ.parameters) : equ.g
_get_ω(equ::LODE) = hasparameters(equ) ? (t,q,v,ω) -> equ.ω(t, q, v, ω, equ.parameters) : equ.ω
_get_v̄(equ::LODE) = hasparameters(equ) ? (t,q,v)   -> equ.v̄(t, q, v, equ.parameters) : equ.v̄
_get_f̄(equ::LODE) = hasparameters(equ) ? (t,q,v,f) -> equ.f̄(t, q, v, f, equ.parameters) : equ.f̄
_get_l(equ::LODE) = hasparameters(equ) ? (t,q,v)   -> equ.lagrangian(t, q, v, equ.parameters) : equ.lagrangian


function get_function_tuple(equ::LODE)
    names = (:ϑ, :f, :g, :l, :ω, :v̄, :f̄)
    equs  = (_get_ϑ(equ), _get_f(equ), _get_g(equ), _get_l(equ), _get_ω(equ), _get_v̄(equ), _get_f̄(equ))

    NamedTuple{names}(equs)
end
