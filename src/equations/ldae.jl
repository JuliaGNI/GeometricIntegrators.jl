@doc raw"""
`LDAE`: Variational Differential Algebraic Equation *EXPERIMENTAL*

Defines an implicit initial value problem
```math
\begin{aligned}
\dot{q} (t) &= v(t) + \lambda(t), &
q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), v(t)) + g(t, q(t), \lambda(t)) + \bar{g} (t, q(t), \mu(t)) , &
p(t_{0}) &= p_{0} , \\
p(t) &= ϑ(t, q(t), v(t)) , \\
0 &= \phi (t, q(t), p(t)) , \\
0 &= \psi (t, q(t), p(t), \dot{q}(t), \dot{p}(t)) ,
\end{aligned}
```
with vector field ``f``, the momentum defined by ``p``, initial conditions ``(q_{0}, p_{0})`` and the solution
``(q,p)`` taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}`` and
the algebraic variables ``(v, \lambda, \mu)`` taking values in
``\mathbb{R}^{d} \times \mathbb{R}^{d} \times \mathbb{R}^{m}``.
This is a special case of a differential algebraic equation with dynamical
variables ``(q,p)`` and algebraic variables ``v``, ``\lambda`` and ``\mu``.

### Parameters

* `DT <: Number`: data type
* `TT <: Real`: time step type
* `AT <: AbstractArray{DT}`: array type
* `ϑType <: Function`: type of `ϑ`
* `fType <: Function`: type of `f`
* `uType <: Function`: type of `u`
* `gType <: Function`: type of `g`
* `ϕType <: Function`: type of `ϕ`
* `ūType <: Function`: type of `ū`
* `ḡType <: Function`: type of `ḡ`
* `ψType <: Function`: type of `ψ`
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
* `u`: function computing the projection for ``q``, for a degenerate system given by ``\lambda``
* `g`: function computing the projection for ``p``, for a degenerate system given by ``\nabla \vartheta (q) \cdot \lambda``
* `ϕ`: primary constraints, for a degenerate system given by ``p - \vartheta (q)``
* `ū`: function computing the secondary projection field ``\bar{u}``, for a degenerate system given by ``\lambda`` (optional)
* `ḡ`: function computing the secondary projection field ``\bar{g}``, for a degenerate system given by ``\lambda \cdot \nabla \vartheta (q)`` (optional)
* `ψ`: secondary constraints, for a degenerate system given by ``\dot{p} - \dot{q} \cdot \nabla \vartheta (q)`` (optional)
* `ω`: function computing the symplectic matrix
* `v̄`: function computing an initial guess for the velocity field ``v`` (optional)
* `f̄`: function computing an initial guess for the force field ``f`` (optional)
* `t₀`: initial time (optional)
* `q₀`: initial condition for dynamical variable `q`
* `p₀`: initial condition for dynamical variable `p`
* `λ₀`: initial condition for algebraic variable `λ`
* `μ₀`: initial condition for algebraic variable `μ` (optional)
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
The funtions `g`, `v̄` and `f̄` are specified by
```julia
    function g(t, q, λ, g)
        g[1] = ...
        g[2] = ...
        ...
    end

    function v̄(t, q, v)
        v[1] = ...
        v[2] = ...
        ...
    end

    function f̄(t, q, v, f)
        f[1] = ...
        f[2] = ...
        ...
    end
```

### Constructors

```julia
LDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, ω, v̄, f̄, t₀, q₀, p₀, λ₀, μ₀, lagrangian, invariants, parameters, periodicity)

LDAE(ϑ, f, u, g, ϕ, l, ω, t₀, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...)
LDAE(ϑ, f, u, g, ϕ, l, ω, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...)
LDAE(ϑ, f, u, g, ϕ, l, ω, t₀, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...)
LDAE(ϑ, f, u, g, ϕ, l, ω, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...)

LDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, l, ω, t₀, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...)
LDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, l, ω, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...)
LDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, l, ω, t₀, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...)
LDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, l, ω, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...)
```

### Keyword arguments

* `v̄ = (t,q,v) -> nothing`
* `f̄ = f`
* `invariants = nothing`
* `parameters = nothing`
* `periodicity = nothing`

"""
struct LDAE{dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
            ϑType <: Function, fType <: Function,
            uType <: Function, gType <: Function, ϕType <: Function,
            ūType <: OptionalFunction, ḡType <: OptionalFunction, ψType <: OptionalFunction,
            ωType <: Function, v̄Type <: Function, f̄Type <: Function,
            lagType <: Function,
            invType <: OptionalNamedTuple,
            parType <: OptionalNamedTuple,
            perType <: OptionalArray{arrayType}} <: AbstractEquationPDAE{dType, tType}

    d::Int
    m::Int

    ϑ::ϑType
    f::fType
    u::uType
    g::gType
    ϕ::ϕType
    ū::ūType
    ḡ::ḡType
    ψ::ψType
    ω::ωType
    v̄::v̄Type
    f̄::f̄Type

    t₀::tType
    q₀::Vector{arrayType}
    p₀::Vector{arrayType}
    λ₀::Vector{arrayType}
    μ₀::Vector{arrayType}

    lagrangian::lagType
    invariants::invType
    parameters::parType
    periodicity::perType

    function LDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, ω, v̄, f̄,
                  t₀::tType,
                  q₀::Vector{arrayType}, p₀::Vector{arrayType},
                  λ₀::Vector{arrayType}, μ₀::Vector{arrayType},
                  lagrangian, invariants, parameters, periodicity) where {
                        dType <: Number, tType <: Real, arrayType <: AbstractArray{dType}}

        d = length(q₀[begin])
        m = length(λ₀[begin])

        @assert 2d ≥ m

        @assert length(q₀) == length(p₀) == length(λ₀) == length(μ₀)

        @assert all(length(q) == d for q in q₀)
        @assert all(length(p) == d for p in p₀)
        @assert all(length(λ) == m for λ in λ₀)
        @assert all(length(μ) == m for μ in μ₀)

        @assert all([ndims(q) == ndims(p) == ndims(λ) == ndims(μ) for (q,p,λ,μ) in zip(q₀,p₀,λ₀,μ₀)])
        
        new{dType, tType, arrayType,
            typeof(ϑ), typeof(f),
            typeof(u), typeof(g), typeof(ϕ),
            typeof(ū), typeof(ḡ), typeof(ψ),
            typeof(ω), typeof(v̄), typeof(f̄),
            typeof(lagrangian), typeof(invariants), typeof(parameters), typeof(periodicity)}(
                d, m, ϑ, f, u, g, ϕ, ū, ḡ, ψ, ω, v̄, f̄, t₀, q₀, p₀, λ₀, μ₀, lagrangian, invariants, parameters, periodicity)
    end
end

_LDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, lagrangian, ω, t₀, q₀, p₀, λ₀, μ₀; v̄=(t,q,v)->nothing, f̄=f, invariants=nothing, parameters=nothing, periodicity=nothing) = LDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, ω, v̄, f̄, t₀, q₀, p₀, λ₀, μ₀, lagrangian, invariants, parameters, periodicity)

LDAE(ϑ, f, u, g, ϕ, l, ω, t₀, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...) = _LDAE(ϑ, f, u, g, ϕ, nothing, nothing, nothing, l, ω, t₀, q₀, p₀, λ₀, μ₀; kwargs...)
LDAE(ϑ, f, u, g, ϕ, l, ω, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...) = LDAE(ϑ, f, u, g, ϕ, l, ω, 0.0, q₀, p₀, λ₀, μ₀; kwargs...)
LDAE(ϑ, f, u, g, ϕ, l, ω, t₀, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...) = LDAE(ϑ, f, u, g, ϕ, l, ω, t₀, [q₀], [p₀], [λ₀], [μ₀]; kwargs...)
LDAE(ϑ, f, u, g, ϕ, l, ω, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...) = LDAE(ϑ, f, u, g, ϕ, l, ω, 0.0, q₀, p₀, λ₀, μ₀; kwargs...)

LDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, l, ω, t₀, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...) = _LDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, l, ω, t₀, q₀, p₀, λ₀, μ₀; kwargs...)
LDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, l, ω, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...) = LDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, l, ω, 0.0, q₀, p₀, λ₀, μ₀; kwargs...)
LDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, l, ω, t₀, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...) = LDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, l, ω, t₀, [q₀], [p₀], [λ₀], [μ₀]; kwargs...)
LDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, l, ω, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...) = LDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, l, ω, 0.0, q₀, p₀, λ₀, μ₀; kwargs...)

const LDAEpsiType{psiT,DT,TT,AT,PT,FT,UT,GT,ΦT,ŪT,ḠT,ΩT,V̄T,F̄T,lagT,invT,parT,perT} = LDAE{DT,TT,AT,PT,FT,UT,GT,ΦT,ŪT,ḠT,psiT,ΩT,V̄T,F̄T,lagT,invT,parT,perT} # type alias for dispatch on secondary constraint type parameter
const LDAEinvType{invT,DT,TT,AT,PT,FT,UT,GT,ΦT,ŪT,ḠT,psiT,ΩT,V̄T,F̄T,lagT,parT,perT} = LDAE{DT,TT,AT,PT,FT,UT,GT,ΦT,ŪT,ḠT,psiT,ΩT,V̄T,F̄T,lagT,invT,parT,perT} # type alias for dispatch on invariants type parameter
const LDAEparType{parT,DT,TT,AT,PT,FT,UT,GT,ΦT,ŪT,ḠT,psiT,ΩT,V̄T,F̄T,lagT,invT,perT} = LDAE{DT,TT,AT,PT,FT,UT,GT,ΦT,ŪT,ḠT,psiT,ΩT,V̄T,F̄T,lagT,invT,parT,perT} # type alias for dispatch on parameters type parameter
const LDAEperType{perT,DT,TT,AT,PT,FT,UT,GT,ΦT,ŪT,ḠT,psiT,ΩT,V̄T,F̄T,lagT,invT,parT} = LDAE{DT,TT,AT,PT,FT,UT,GT,ΦT,ŪT,ḠT,psiT,ΩT,V̄T,F̄T,lagT,invT,parT,perT} # type alias for dispatch on periodicity type parameter

Base.hash(dae::LDAE, h::UInt) = hash(dae.d, hash(dae.m,
          hash(dae.ϑ, hash(dae.f, 
          hash(dae.u, hash(dae.g, hash(dae.ϕ, 
          hash(dae.ū, hash(dae.ḡ, hash(dae.ψ,
          hash(dae.ω, hash(dae.v̄, hash(dae.f̄,
          hash(dae.t₀, hash(dae.q₀, hash(dae.p₀, hash(dae.λ₀, hash(dae.μ₀,
          hash(dae.invariants, hash(dae.parameters, hash(dae.periodicity, h)))))))))))))))))))))

Base.:(==)(dae1::LDAE, dae2::LDAE) = (
                                dae1.d == dae2.d
                             && dae1.m == dae2.m
                             && dae1.ϑ == dae2.ϑ
                             && dae1.f == dae2.f
                             && dae1.u == dae2.u
                             && dae1.g == dae2.g
                             && dae1.ϕ == dae2.ϕ
                             && dae1.ū == dae2.ū
                             && dae1.ḡ == dae2.ḡ
                             && dae1.ψ == dae2.ψ
                             && dae1.ω == dae2.ω
                             && dae1.v̄ == dae2.v̄
                             && dae1.f̄ == dae2.f̄
                             && dae1.t₀ == dae2.t₀
                             && dae1.q₀ == dae2.q₀
                             && dae1.p₀ == dae2.p₀
                             && dae1.λ₀ == dae2.λ₀
                             && dae1.μ₀ == dae2.μ₀
                             && dae1.lagrangian == dae2.lagrangian
                             && dae1.invariants  == dae2.invariants
                             && dae1.parameters  == dae2.parameters
                             && dae1.periodicity == dae2.periodicity)

function Base.similar(equ::LDAE, t₀::Real, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector; parameters=equ.parameters)
    @assert all([length(q) == ndims(equ) for q in q₀])
    @assert all([length(p) == ndims(equ) for p in p₀])
    @assert all([length(λ) == ndims(equ) for λ in λ₀])
    @assert all([length(μ) == ndims(equ) for μ in μ₀])
    LDAE(equ.ϑ, equ.f, equ.u, equ.g, equ.ϕ, equ.ū, equ.ḡ, equ.ψ, equ.lagrangian, equ.ω, t₀, q₀, p₀, λ₀, μ₀;
         v̄=equ.v̄, f̄=equ.f̄, invariants=equ.invariants, parameters=parameters, periodicity=equ.periodicity)
end

Base.similar(equ::LDAE, q₀, p₀, λ₀=get_λ₀(q₀, equ.λ₀), μ₀=get_λ₀(λ₀, equ.μ₀); kwargs...) = similar(equ, equ.t₀, q₀, p₀, λ₀, μ₀; kwargs...)
Base.similar(equ::LDAE, t₀::Real, q₀::State, p₀::State, λ₀::State=get_λ₀(q₀, equ.λ₀), μ₀=get_λ₀(λ₀, equ.μ₀); kwargs...) = similar(equ, t₀, [q₀], [p₀], [λ₀], [μ₀]; kwargs...)

hassecondary(::LDAEpsiType{<:Nothing}) = false
hassecondary(::LDAEpsiType{<:Function}) = true

hasinvariants(::LDAEinvType{<:Nothing}) = false
hasinvariants(::LDAEinvType{<:NamedTuple}) = true

hasparameters(::LDAEparType{<:Nothing}) = false
hasparameters(::LDAEparType{<:NamedTuple}) = true

hasperiodicity(::LDAEperType{<:Nothing}) = false
hasperiodicity(::LDAEperType{<:AbstractArray}) = true

@inline Base.axes(equation::LDAE) = axes(equation.q₀[begin])
@inline Base.ndims(equation::LDAE) = equation.d
@inline Common.nsamples(equ::LDAE) = length(eachindex(equ.q₀))
@inline Common.nconstraints(equation::LDAE) = equation.m

@inline Common.periodicity(equation::LDAE) = hasperiodicity(equation) ? equation.periodicity : zero(equation.q₀[begin])
@inline initial_conditions(equation::LDAE) = (equation.t₀, equation.q₀, equation.p₀, equation.λ₀, equation.μ₀)

_get_ϑ(equ::LDAE) = hasparameters(equ) ? (t,q,v,ϑ)     -> equ.ϑ(t, q, v, ϑ, equ.parameters) : equ.ϑ
_get_f(equ::LDAE) = hasparameters(equ) ? (t,q,v,f)     -> equ.f(t, q, v, f, equ.parameters) : equ.f
_get_u(equ::LDAE) = hasparameters(equ) ? (t,q,v,u)     -> equ.u(t, q, v, u, equ.parameters) : equ.u
_get_g(equ::LDAE) = hasparameters(equ) ? (t,q,v,g)     -> equ.g(t, q, v, g, equ.parameters) : equ.g
_get_ϕ(equ::LDAE) = hasparameters(equ) ? (t,q,v,ϕ)     -> equ.ϕ(t, q, v, ϕ, equ.parameters) : equ.ϕ
_get_ū(equ::LDAE) = hasparameters(equ) ? (t,q,λ,ū)     -> equ.ū(t, q, λ, ū, equ.parameters) : equ.ū
_get_ḡ(equ::LDAE) = hasparameters(equ) ? (t,q,λ,ḡ)     -> equ.ḡ(t, q, λ, ḡ, equ.parameters) : equ.ḡ
_get_ψ(equ::LDAE) = hasparameters(equ) ? (t,q,p,v,f,ψ) -> equ.ψ(t, q, p, v, f, ψ, equ.parameters) : equ.ψ
_get_ω(equ::LDAE) = hasparameters(equ) ? (t,q,v,ω)     -> equ.ω(t, q, v, ω, equ.parameters) : equ.ω
_get_v̄(equ::LDAE) = hasparameters(equ) ? (t,q,v)       -> equ.v̄(t, q, v, equ.parameters) : equ.v̄
_get_f̄(equ::LDAE) = hasparameters(equ) ? (t,q,v,f)     -> equ.f̄(t, q, v, f, equ.parameters) : equ.f̄
_get_l(equ::LDAE) = hasparameters(equ) ? (t,q,v)       -> equ.lagrangian(t, q, v, equ.parameters) : equ.lagrangian


function get_function_tuple(equ::LDAE)
    names = (:ϑ, :f, :u, :g, :ϕ)
    equs  = (_get_ϑ(equ), _get_f(equ), _get_u(equ), _get_g(equ), _get_ϕ(equ))

    if hassecondary(equ)
        names = (names..., :ū, :ḡ, :ψ)
        equs  = (equs..., _get_ū(equ), _get_ḡ(equ), _get_ψ(equ))
    end

    names = (names..., :l, :ω, :v̄, :f̄)
    equs  = (equs..., _get_l(equ), _get_ω(equ), _get_v̄(equ), _get_f̄(equ))

    NamedTuple{names}(equs)
end
