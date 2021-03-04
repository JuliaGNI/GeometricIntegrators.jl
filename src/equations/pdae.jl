@doc raw"""
`PDAE`: Partitioned Differential Algebraic Equation

Defines a partitioned differential algebraic initial value problem
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t), p(t)) + u(t, q(t), p(t), \lambda(t)) , & q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), p(t)) + r(t, q(t), p(t), \lambda(t)) , & p(t_{0}) &= p_{0} , \\
0 &= \phi (t, q(t), p(t), \lambda(t)) , & \lambda(t_{0}) &= \lambda_{0} ,
\end{aligned}
```
with vector fields ``v`` and ``f``, projection ``u`` and ``r``,
algebraic constraint ``\phi=0``,
conditions ``(q_{0}, p_{0})`` and ``\lambda_{0}``, the dynamical variables
``(q,p)`` taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}`` and
the algebraic variable ``\lambda`` taking values in ``\mathbb{R}^{n}``.

### Parameters

* `DT <: Number`: data type
* `TT <: Real`: time step type
* `AT <: AbstractArray{DT}`: array type
* `vType <: Function`: type of `v`
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
* `v`: function computing the vector field ``v``
* `f`: function computing the vector field ``f``
* `u`: function computing the projection for ``q``
* `g`: function computing the projection for ``p``
* `ϕ`: algebraic constraints
* `ū`: function computing the secondary projection field ``\bar{u}`` (optional)
* `ḡ`: function computing the secondary projection field ``\bar{g}`` (optional)
* `ψ`: secondary constraints (optional)
* `v̄`: function computing an initial guess for the velocity field ``v`` (optional)
* `f̄`: function computing an initial guess for the force field ``f`` (optional)
* `t₀`: initial time
* `q₀`: initial condition for dynamical variable ``q``
* `p₀`: initial condition for dynamical variable ``p``
* `λ₀`: initial condition for algebraic variable ``\lambda``
* `μ₀`: initial condition for algebraic variable ``μ`` (optional)
* `invariants`: either a `NamedTuple` containing the equation's invariants or `nothing`
* `parameters`: either a `NamedTuple` containing the equation's parameters or `nothing`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions

### Constructors

```julia
PDAE(v, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, t₀, q₀, p₀, λ₀, invariants, parameters, periodicity)

PDAE(v, f, u, g, ϕ, t₀, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...)
PDAE(v, f, u, g, ϕ, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...)
PDAE(v, f, u, g, ϕ, t₀, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...)
PDAE(v, f, u, g, ϕ, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...)

PDAE(v, f, u, g, ϕ, ū, ḡ, ψ, t₀, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...)
PDAE(v, f, u, g, ϕ, ū, ḡ, ψ, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...)
PDAE(v, f, u, g, ϕ, ū, ḡ, ψ, t₀, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...)
PDAE(v, f, u, g, ϕ, ū, ḡ, ψ, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...)
```

"""
struct PDAE{dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
            vType <: Function, fType <: Function,
            uType <: Function, gType <: Function, ϕType <: Function,
            ūType <: OptionalFunction, ḡType <: OptionalFunction, ψType <: OptionalFunction,
            v̄Type <: Function, f̄Type <: Function,
            invType <: OptionalNamedTuple,
            parType <: OptionalNamedTuple,
            perType <: OptionalArray{arrayType}} <: AbstractEquationPDAE{dType, tType}

    d::Int
    m::Int

    v::vType
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

    function PDAE(v::vType, f::fType, u::uType, g::gType, ϕ::ϕType, ū::ūType, ḡ::ḡType, ψ::ψType, v̄::v̄Type, f̄::f̄Type,
                  t₀::tType, q₀::Vector{arrayType}, p₀::Vector{arrayType}, λ₀::Vector{arrayType}, μ₀::Vector{arrayType},
                  invariants::invType, parameters::parType, periodicity::perType) where {
                        dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
                        vType <: Function, fType <: Function,
                        uType <: Function, gType <: Function, ϕType <: Function,
                        ūType <: OptionalFunction, ḡType <: OptionalFunction, ψType <: OptionalFunction,
                        v̄Type <: Function, f̄Type <: Function, 
                        invType <: OptionalNamedTuple,
                        parType <: OptionalNamedTuple,
                        perType <: OptionalArray{arrayType}}

        d = length(q₀[begin])
        m = length(λ₀[begin])

        @assert 2d ≥ m

        @assert length(q₀) == length(p₀) == length(λ₀)

        @assert all(length(q) == d for q in q₀)
        @assert all(length(p) == d for p in p₀)
        @assert all(length(λ) == m for λ in λ₀)
        @assert all(length(μ) == m for μ in μ₀)

        @assert all([ndims(q) == ndims(p) == ndims(λ) == ndims(μ) for (q,p,λ,μ) in zip(q₀,p₀,λ₀,μ₀)])

        new{dType, tType, arrayType, vType, fType, uType, gType, ϕType, ūType, ḡType, ψType, v̄Type, f̄Type, invType, parType, perType}(d, m, v, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, t₀, q₀, p₀, λ₀, μ₀, invariants, parameters, periodicity)
    end
end

_PDAE(v, f, u, g, ϕ, ū, ḡ, ψ, t₀, q₀, p₀, λ₀, μ₀; v̄=v, f̄=f, invariants=nothing, parameters=nothing, periodicity=nothing) = PDAE(v, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, t₀, q₀, p₀, λ₀, μ₀, invariants, parameters, periodicity)

PDAE(v, f, u, g, ϕ, t₀, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...) = _PDAE(v, f, u, g, ϕ, nothing, nothing, nothing, t₀, q₀, p₀, λ₀, μ₀; kwargs...)
PDAE(v, f, u, g, ϕ, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...) = PDAE(v, f, u, g, ϕ, 0.0, q₀, p₀, λ₀, μ₀; kwargs...)
PDAE(v, f, u, g, ϕ, t₀, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...) = PDAE(v, f, u, g, ϕ, t₀, [q₀], [p₀], [λ₀], [μ₀]; kwargs...)
PDAE(v, f, u, g, ϕ, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...) = PDAE(v, f, u, g, ϕ, 0.0, q₀, p₀, λ₀, μ₀; kwargs...)

PDAE(v, f, u, g, ϕ, ū, ḡ, ψ, t₀, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...) = _PDAE(v, f, u, g, ϕ, ū, ḡ, ψ, t₀, q₀, p₀, λ₀, μ₀; kwargs...)
PDAE(v, f, u, g, ϕ, ū, ḡ, ψ, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...) = PDAE(v, f, u, g, ϕ, ū, ḡ, ψ, 0.0, q₀, p₀, λ₀, μ₀; kwargs...)
PDAE(v, f, u, g, ϕ, ū, ḡ, ψ, t₀, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...) = PDAE(v, f, u, g, ϕ, ū, ḡ, ψ, t₀, [q₀], [p₀], [λ₀], [μ₀]; kwargs...)
PDAE(v, f, u, g, ϕ, ū, ḡ, ψ, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...) = PDAE(v, f, u, g, ϕ, ū, ḡ, ψ, 0.0, q₀, p₀, λ₀, μ₀; kwargs...)

const PDAEpsiType{psiT,DT,TT,AT,VT,FT,UT,GT,ΦT,ŪT,ḠT,V̄T,F̄T,invT,parT,perT} = PDAE{DT,TT,AT,VT,FT,UT,GT,ΦT,ŪT,ḠT,psiT,V̄T,F̄T,invT,parT,perT} # type alias for dispatch on secondary constraint type parameter
const PDAEinvType{invT,DT,TT,AT,VT,FT,UT,GT,ΦT,ŪT,ḠT,psiT,V̄T,F̄T,parT,perT} = PDAE{DT,TT,AT,VT,FT,UT,GT,ΦT,ŪT,ḠT,psiT,V̄T,F̄T,invT,parT,perT} # type alias for dispatch on invariants type parameter
const PDAEparType{parT,DT,TT,AT,VT,FT,UT,GT,ΦT,ŪT,ḠT,psiT,V̄T,F̄T,invT,perT} = PDAE{DT,TT,AT,VT,FT,UT,GT,ΦT,ŪT,ḠT,psiT,V̄T,F̄T,invT,parT,perT} # type alias for dispatch on parameters type parameter
const PDAEperType{perT,DT,TT,AT,VT,FT,UT,GT,ΦT,ŪT,ḠT,psiT,V̄T,F̄T,invT,parT} = PDAE{DT,TT,AT,VT,FT,UT,GT,ΦT,ŪT,ḠT,psiT,V̄T,F̄T,invT,parT,perT} # type alias for dispatch on periodicity type parameter

Base.hash(dae::PDAE, h::UInt) = hash(dae.d, hash(dae.m,
        hash(dae.v, hash(dae.f, 
        hash(dae.u, hash(dae.g, hash(dae.ϕ,
        hash(dae.ū, hash(dae.ḡ, hash(dae.ψ,
        hash(dae.v̄, hash(dae.f̄,
        hash(dae.t₀, hash(dae.q₀, hash(dae.p₀, hash(dae.λ₀, hash(dae.μ₀,
        hash(dae.invariants, hash(dae.parameters, hash(dae.periodicity, h))))))))))))))))))))

Base.:(==)(dae1::PDAE, dae2::PDAE) = (
                                dae1.d == dae2.d
                             && dae1.m == dae2.m
                             && dae1.v == dae2.v
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

function Base.similar(equ::PDAE, t₀::Real, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector; parameters=equ.parameters)
    @assert all([length(q) == equ.d for q in q₀])
    @assert all([length(p) == equ.d for p in p₀])
    @assert all([length(λ) == equ.m for λ in λ₀])
    @assert all([length(μ) == equ.m for μ in μ₀])
    PDAE(equ.v, equ.f, equ.u, equ.g, equ.ϕ, equ.ū, equ.ḡ, equ.ψ, t₀, q₀, p₀, λ₀, μ₀; v̄=equ.v̄, f̄=equ.f̄,
         invariants=equ.invariants, parameters=parameters, periodicity=equ.periodicity)
end

Base.similar(equ::PDAE, q₀, p₀, λ₀=get_λ₀(q₀, equ.λ₀), μ₀=get_λ₀(λ₀, equ.μ₀); kwargs...) = similar(equ, equ.t₀, q₀, p₀, λ₀, μ₀; kwargs...)
Base.similar(equ::PDAE, t₀::Real, q₀::State, p₀::State, λ₀::State=get_λ₀(q₀, equ.λ₀), μ₀::State=get_λ₀(λ₀, equ.μ₀); kwargs...) = similar(equ, t₀, [q₀], [p₀], [λ₀], [μ₀]; kwargs...)

hassecondary(::PDAEpsiType{<:Nothing}) = false
hassecondary(::PDAEpsiType{<:Function}) = true

hasinvariants(::PDAEinvType{<:Nothing}) = false
hasinvariants(::PDAEinvType{<:NamedTuple}) = true

hasparameters(::PDAEparType{<:Nothing}) = false
hasparameters(::PDAEparType{<:NamedTuple}) = true

hasperiodicity(::PDAEperType{<:Nothing}) = false
hasperiodicity(::PDAEperType{<:AbstractArray}) = true

@inline Base.axes(equation::PDAE) = axes(equation.q₀[begin])
@inline Base.ndims(equation::PDAE) = equation.d
@inline Common.nsamples(equation::PDAE) = length(eachindex(equation.q₀))
@inline Common.nconstraints(equation::PDAE) = equation.m

@inline Common.periodicity(equation::PDAE) = hasperiodicity(equation) ? equation.periodicity : zero(equation.q₀[begin])
@inline initial_conditions(equation::PDAE) = (equation.t₀, equation.q₀, equation.p₀, equation.λ₀, equation.μ₀)

_get_v(equ::PDAE) = hasparameters(equ) ? (t,q,p,v)     -> equ.v(t, q, p, v, equ.parameters) : equ.v
_get_f(equ::PDAE) = hasparameters(equ) ? (t,q,p,f)     -> equ.f(t, q, p, f, equ.parameters) : equ.f
_get_u(equ::PDAE) = hasparameters(equ) ? (t,q,p,λ,u)   -> equ.u(t, q, p, λ, u, equ.parameters) : equ.u
_get_g(equ::PDAE) = hasparameters(equ) ? (t,q,p,λ,g)   -> equ.g(t, q, p, λ, g, equ.parameters) : equ.g
_get_ϕ(equ::PDAE) = hasparameters(equ) ? (t,q,p,ϕ)     -> equ.ϕ(t, q, p, ϕ, equ.parameters) : equ.ϕ
_get_ū(equ::PDAE) = hasparameters(equ) ? (t,q,p,λ,u)   -> equ.ū(t, q, p, λ, u, equ.parameters) : equ.ū
_get_ḡ(equ::PDAE) = hasparameters(equ) ? (t,q,p,λ,g)   -> equ.ḡ(t, q, p, λ, g, equ.parameters) : equ.ḡ
_get_ψ(equ::PDAE) = hasparameters(equ) ? (t,q,p,v,f,ψ) -> equ.ψ(t, q, p, v, f, ψ, equ.parameters) : equ.ψ
_get_v̄(equ::PDAE) = hasparameters(equ) ? (t,q,p,v)     -> equ.v̄(t, q, p, v, equ.parameters) : equ.v̄
_get_f̄(equ::PDAE) = hasparameters(equ) ? (t,q,p,f)     -> equ.f̄(t, q, p, f, equ.parameters) : equ.f̄


function get_function_tuple(equ::PDAE)
    names = (:v,:f,:u,:g,:ϕ)
    equs  = (_get_v(equ), _get_f(equ), _get_u(equ), _get_g(equ), _get_ϕ(equ))

    if hassecondary(equ)
        names = (names..., :ū, :ḡ, :ψ)
        equs  = (equs..., _get_ū(equ), _get_ḡ(equ), _get_ψ(equ))
    end

    names = (names..., :v̄, :f̄)
    equs  = (equs..., _get_v̄(equ), _get_f̄(equ))
    
    NamedTuple{names}(equs)
end
