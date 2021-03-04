@doc raw"""
`IDAE`: Implicit Differential Algebraic Equation

Defines a partitioned differential algebraic initial value problem
```math
\begin{aligned}
\dot{q} (t) &= v(t) + u(t, q(t), p(t), \lambda(t)) , & q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), v(t)) + r(t, q(t), p(t), \lambda(t)) , & p(t_{0}) &= p_{0} , \\
p(t) &= p(t, q(t), v(t)) , && \\
0 &= \phi (t, q(t), p(t), \lambda(t)) , & \lambda(t_{0}) &= \lambda_{0} ,
\end{aligned}
```
with vector field ``f``, the momentum defined by ``p``, projection ``u`` and ``r``,
algebraic constraint ``\phi=0``,
conditions ``(q_{0}, p_{0})`` and ``\lambda_{0}``, the dynamical variables
``(q,p)`` taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}`` and
the algebraic variable ``\lambda`` taking values in ``\mathbb{R}^{n}``.

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
* `v̄Type <: Function`: type of `v̄`
* `f̄Type <: Function`: type of `f̄`
* `invType <: OptionalNamedTuple`: invariants type
* `parType <: OptionalNamedTuple`: parameters type
* `perType <: OptionalArray{AT}`: periodicity type

### Fields

* `d`: dimension of dynamical variables ``q`` and ``p`` as well as the vector fields ``f`` and ``p``
* `m`: dimension of algebraic variable ``\lambda`` and the constraint ``\phi``
* `ϑ`: function determining the momentum
* `f`: function computing the vector field ``f``
* `u`: function computing the projection for ``q``
* `g`: function computing the projection for ``p``
* `ϕ`: algebraic constraints
* `ū`: function computing the secondary projection field ``\bar{u}`` (optional)
* `ḡ`: function computing the secondary projection field ``\bar{g}`` (optional)
* `ψ`: secondary constraints (optional)
* `v̄`: function computing an initial guess for the velocity field ``v`` (optional)
* `f̄`: function computing an initial guess for the force field ``f`` (optional)
* `t₀`: initial time (optional)
* `q₀`: initial condition for dynamical variable ``q``
* `p₀`: initial condition for dynamical variable ``p``
* `λ₀`: initial condition for algebraic variable `λ`
* `μ₀`: initial condition for algebraic variable `μ` (optional)
* `invariants`: either a `NamedTuple` containing the equation's invariants or `nothing`
* `parameters`: either a `NamedTuple` containing the equation's parameters or `nothing`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions

### Constructors

```julia
IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, t₀, q₀, p₀, λ₀, μ₀, invariants, parameters, periodicity)

IDAE(ϑ, f, u, g, ϕ, t₀, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...)
IDAE(ϑ, f, u, g, ϕ, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...)
IDAE(ϑ, f, u, g, ϕ, t₀, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...)
IDAE(ϑ, f, u, g, ϕ, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...)

IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, t₀, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...)
IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...)
IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, t₀, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...)
IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...)
```

### Keyword arguments

* `v̄ = (t,q,v) -> nothing`
* `f̄ = f`
* `invariants = nothing`
* `parameters = nothing`
* `periodicity = nothing`


"""
struct IDAE{dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
            ϑType <: Function, fType <: Function,
            uType <: Function, gType <: Function, ϕType <: Function,
            ūType <: OptionalFunction, ḡType <: OptionalFunction, ψType <: OptionalFunction,
            v̄Type <: Function, f̄Type <: Function,
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
    v̄::v̄Type
    f̄::f̄Type
    t₀::tType
    q₀::Vector{arrayType}
    p₀::Vector{arrayType}
    λ₀::Vector{arrayType}
    μ₀::Vector{arrayType}
    invariants::invType
    parameters::parType
    periodicity::perType

    function IDAE(ϑ::ϑType, f::fType, u::uType, g::gType, ϕ::ϕType, ū::ūType, ḡ::ḡType, ψ::ψType, v̄::v̄Type, f̄::f̄Type,
                  t₀::tType, q₀::Vector{arrayType}, p₀::Vector{arrayType}, λ₀::Vector{arrayType}, μ₀::Vector{arrayType},
                  invariants::invType, parameters::parType, periodicity::perType) where {
                        dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
                        ϑType <: Function, fType <: Function,
                        uType <: Function, gType <: Function, ϕType <: Function,
                        ūType <: OptionalFunction, ḡType <: OptionalFunction, ψType <: OptionalFunction,
                        v̄Type <: Function, f̄Type <: Function,
                        invType <: OptionalNamedTuple,
                        parType <: OptionalNamedTuple,
                        perType <: OptionalArray{arrayType}}

        d = length(q₀[begin])
        m = length(λ₀[begin])

        @assert 2d ≥ m

        @assert length(q₀) == length(p₀) == length(λ₀) == length(μ₀)

        @assert all(length(q) == d for q in q₀)
        @assert all(length(p) == d for p in p₀)
        @assert all(length(λ) == m for λ in λ₀)
        @assert all(length(μ) == m for μ in μ₀)

        @assert all([ndims(q) == ndims(p) == ndims(λ) == ndims(μ) for (q,p,λ,μ) in zip(q₀,p₀,λ₀,μ₀)])

        new{dType, tType, arrayType, ϑType, fType, uType, gType, ϕType, ūType, ḡType, ψType, v̄Type, f̄Type, invType, parType, perType}(d, m, ϑ, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, t₀, q₀, p₀, λ₀, μ₀, invariants, parameters, periodicity)
    end
end

_IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, t₀, q₀, p₀, λ₀, μ₀; v̄=(t,q,v)->nothing, f̄=f, invariants=nothing, parameters=nothing, periodicity=nothing) = IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, t₀, q₀, p₀, λ₀, μ₀, invariants, parameters, periodicity)

IDAE(ϑ, f, u, g, ϕ, t₀, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...) = _IDAE(ϑ, f, u, g, ϕ, nothing, nothing, nothing, t₀, q₀, p₀, λ₀, μ₀; kwargs...)
IDAE(ϑ, f, u, g, ϕ, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...) = IDAE(ϑ, f, u, g, ϕ, 0.0, q₀, p₀, λ₀, μ₀; kwargs...)
IDAE(ϑ, f, u, g, ϕ, t₀, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...) = IDAE(ϑ, f, u, g, ϕ, t₀, [q₀], [p₀], [λ₀], [μ₀]; kwargs...)
IDAE(ϑ, f, u, g, ϕ, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...) = IDAE(ϑ, f, u, g, ϕ, 0.0, q₀, p₀, λ₀, μ₀; kwargs...)

IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, t₀, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...) = _IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, t₀, q₀, p₀, λ₀, μ₀; kwargs...)
IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...) = IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, 0.0, q₀, p₀, λ₀, μ₀; kwargs...)
IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, t₀, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...) = IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, t₀, [q₀], [p₀], [λ₀], [μ₀]; kwargs...)
IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...) = IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, 0.0, q₀, p₀, λ₀, μ₀; kwargs...)

const IDAEpsiType{psiT,DT,TT,AT,PT,FT,UT,GT,ΦT,ŪT,ḠT,V̄T,F̄T,invT,parT,perT} = IDAE{DT,TT,AT,PT,FT,UT,GT,ΦT,ŪT,ḠT,psiT,V̄T,F̄T,invT,parT,perT} # type alias for dispatch on secondary constraint type parameter
const IDAEinvType{invT,DT,TT,AT,PT,FT,UT,GT,ΦT,ŪT,ḠT,psiT,V̄T,F̄T,parT,perT} = IDAE{DT,TT,AT,PT,FT,UT,GT,ΦT,ŪT,ḠT,psiT,V̄T,F̄T,invT,parT,perT} # type alias for dispatch on invariants type parameter
const IDAEparType{parT,DT,TT,AT,PT,FT,UT,GT,ΦT,ŪT,ḠT,psiT,V̄T,F̄T,invT,perT} = IDAE{DT,TT,AT,PT,FT,UT,GT,ΦT,ŪT,ḠT,psiT,V̄T,F̄T,invT,parT,perT} # type alias for dispatch on parameters type parameter
const IDAEperType{perT,DT,TT,AT,PT,FT,UT,GT,ΦT,ŪT,ḠT,psiT,V̄T,F̄T,invT,parT} = IDAE{DT,TT,AT,PT,FT,UT,GT,ΦT,ŪT,ḠT,psiT,V̄T,F̄T,invT,parT,perT} # type alias for dispatch on periodicity type parameter

Base.hash(dae::IDAE, h::UInt) = hash(dae.d, hash(dae.m,
        hash(dae.ϑ, hash(dae.f, 
        hash(dae.u, hash(dae.g, hash(dae.ϕ, 
        hash(dae.ū, hash(dae.ḡ, hash(dae.ψ,
        hash(dae.v̄, hash(dae.f̄,
        hash(dae.t₀, hash(dae.q₀, hash(dae.p₀, hash(dae.λ₀, hash(dae.μ₀,
        hash(dae.invariants, hash(dae.parameters, hash(dae.periodicity, h))))))))))))))))))))

Base.:(==)(dae1::IDAE, dae2::IDAE) = (
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
                             && dae1.v̄ == dae2.v̄
                             && dae1.f̄ == dae2.f̄
                             && dae1.t₀ == dae2.t₀
                             && dae1.q₀ == dae2.q₀
                             && dae1.p₀ == dae2.p₀
                             && dae1.λ₀ == dae2.λ₀
                             && dae1.μ₀ == dae2.μ₀
                             && dae1.invariants  == dae2.invariants
                             && dae1.parameters  == dae2.parameters
                             && dae1.periodicity == dae2.periodicity)

function Base.similar(equ::IDAE, t₀::Real, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector; parameters=equ.parameters)
    @assert all([length(q) == ndims(equ) for q in q₀])
    @assert all([length(p) == ndims(equ) for p in p₀])
    @assert all([length(λ) == equ.m for λ in λ₀])
    IDAE(equ.ϑ, equ.f, equ.u, equ.g, equ.ϕ, equ.ū, equ.ḡ, equ.ψ, t₀, q₀, p₀, λ₀, μ₀; v̄=equ.v̄, f̄=equ.f̄,
         invariants=equ.invariants, parameters=parameters, periodicity=equ.periodicity)
end

Base.similar(equ::IDAE, q₀, p₀, λ₀=get_λ₀(q₀, equ.λ₀), μ₀=get_λ₀(λ₀, equ.μ₀); kwargs...) = similar(equ, equ.t₀, q₀, p₀, λ₀, μ₀; kwargs...)
Base.similar(equ::IDAE, t₀::Real, q₀::State, p₀::State, λ₀::State=get_λ₀(q₀, equ.λ₀), μ₀=get_λ₀(λ₀, equ.μ₀); kwargs...) = similar(equ, t₀, [q₀], [p₀], [λ₀], [μ₀]; kwargs...)

hassecondary(::IDAEpsiType{<:Nothing}) = false
hassecondary(::IDAEpsiType{<:Function}) = true

hasinvariants(::IDAEinvType{<:Nothing}) = false
hasinvariants(::IDAEinvType{<:NamedTuple}) = true

hasparameters(::IDAEparType{<:Nothing}) = false
hasparameters(::IDAEparType{<:NamedTuple}) = true

hasperiodicity(::IDAEperType{<:Nothing}) = false
hasperiodicity(::IDAEperType{<:AbstractArray}) = true

@inline Base.axes(equation::IDAE) = axes(equation.q₀[begin])
@inline Base.ndims(equation::IDAE) = equation.d
@inline Common.nsamples(equation::IDAE) = length(eachindex(equation.q₀))
@inline Common.nconstraints(equation::IDAE) = equation.m

@inline Common.periodicity(equation::IDAE) = hasperiodicity(equation) ? equation.periodicity : zero(equation.q₀[begin])
@inline initial_conditions(equation::IDAE) = (equation.t₀, equation.q₀, equation.p₀, equation.λ₀, equation.μ₀)

_get_ϑ(equ::IDAE) = hasparameters(equ) ? (t,q,v,ϑ)     -> equ.ϑ(t, q, v, ϑ, equ.parameters) : equ.ϑ
_get_f(equ::IDAE) = hasparameters(equ) ? (t,q,v,f)     -> equ.f(t, q, v, f, equ.parameters) : equ.f
_get_u(equ::IDAE) = hasparameters(equ) ? (t,q,p,λ,g)   -> equ.u(t, q, p, λ, g, equ.parameters) : equ.u
_get_g(equ::IDAE) = hasparameters(equ) ? (t,q,p,λ,g)   -> equ.g(t, q, p, λ, g, equ.parameters) : equ.g
_get_ϕ(equ::IDAE) = hasparameters(equ) ? (t,q,p,ϕ)     -> equ.ϕ(t, q, p, ϕ, equ.parameters) : equ.ϕ
_get_ū(equ::IDAE) = hasparameters(equ) ? (t,q,p,μ,g)   -> equ.ū(t, q, p, μ, g, equ.parameters) : equ.ū
_get_ḡ(equ::IDAE) = hasparameters(equ) ? (t,q,p,μ,g)   -> equ.ḡ(t, q, p, μ, g, equ.parameters) : equ.ḡ
_get_ψ(equ::IDAE) = hasparameters(equ) ? (t,q,p,v,f,ψ) -> equ.ψ(t, q, p, v, f, ψ, equ.parameters) : equ.ψ
_get_v̄(equ::IDAE) = hasparameters(equ) ? (t,q,v)       -> equ.v̄(t, q, v, equ.parameters) : equ.v̄
_get_f̄(equ::IDAE) = hasparameters(equ) ? (t,q,v,f)     -> equ.f̄(t, q, v, f, equ.parameters) : equ.f̄
_get_h(equ::IDAE) = hasparameters(equ) ? (t,q)         -> equ.h(t, q, equ.parameters) : equ.h

function get_function_tuple(equ::IDAE)
    names = (:ϑ,:f,:u,:g,:ϕ)
    equs  = (_get_ϑ(equ), _get_f(equ), _get_u(equ), _get_g(equ), _get_ϕ(equ))

    if hassecondary(equ)
        names = (names..., :ū, :ḡ, :ψ)
        equs  = (equs..., _get_ū(equ), _get_ḡ(equ), _get_ψ(equ))
    end

    names = (names..., :v̄, :f̄)
    equs  = (equs..., _get_v̄(equ), _get_f̄(equ))

    NamedTuple{names}(equs)
end
