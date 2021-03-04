@doc raw"""
`DAE`: Differential Algebraic Equation

Defines a differential algebraic initial value problem
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t)) + u(t, q(t), \lambda(t)) , & q(t_{0}) &= q_{0} , \\
0 &= \phi (t, q(t), \lambda(t)) , & \lambda(t_{0}) &= \lambda_{0} ,
\end{aligned}
```
with vector field ``v``, projection ``u``, algebraic constraint ``\phi=0``,
initial conditions ``q_{0}`` and ``\lambda_{0}``, the dynamical variable ``q``
taking values in ``\mathbb{R}^{m}`` and the algebraic variable ``\lambda``
taking values in ``\mathbb{R}^{n}``.

### Parameters

* `DT <: Number`: data type
* `TT <: Real`: time step type
* `AT <: AbstractArray{DT}`: array type
* `vType <: Function`: type of `v`
* `uType <: Function`: type of `u`
* `ūType <: OptionalFunction`: type of `ū`
* `ϕType <: Function`: type of `ϕ`
* `ψType <: OptionalFunction`: type of `ψ`
* `v̄Type <: Function`: type of `v̄`
* `invType <: OptionalNamedTuple`: invariants type
* `parType <: OptionalNamedTuple`: parameters type
* `perType <: OptionalArray{AT}`: periodicity type

### Fields

* `d`: dimension of dynamical variable ``q`` and the vector field ``v``
* `m`: dimension of algebraic variable ``\lambda`` and the constraint ``\phi``
* `v`: function computing the vector field
* `u`: function computing the projection
* `ū`: function computing the secondary projection field ``\bar{u}`` (optional)
* `ϕ`: algebraic constraint
* `ψ`: secondary constraints (optional)
* `v̄`: function computing an initial guess for the velocity field ``v`` (defaults to `v`)
* `t₀`: initial time
* `q₀`: initial condition for dynamical variable ``q``
* `λ₀`: initial condition for algebraic variable ``\lambda``
* `invariants`: either a `NamedTuple` containing the equation's invariants or `nothing`
* `parameters`: either a `NamedTuple` containing the equation's parameters or `nothing`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions

The function `v`, providing the vector field, takes three arguments,
`v(t, q, v)`, the functions `u` and `ϕ`, providing the projection and the
algebraic constraint take four arguments, `u(t, q, λ, u)` and `ϕ(t, q, λ, ϕ)`,
where `t` is the current time, `q` and `λ` are the current solution vectors,
and `v`, `u` and `ϕ` are the vectors which hold the result of evaluating the
vector field ``v``, the projection ``u`` and the algebraic constraint ``\phi``
on `t`, `q` and `λ`.

### Constructors

```julia
DAE(v, u, ū, ϕ, ψ, v̄, t₀, q₀, λ₀, invariants, parameters, periodicity)

DAE(v, u, ϕ, t₀, q₀::StateVector, λ₀::StateVector; kwargs...)
DAE(v, u, ϕ, q₀::StateVector, λ₀::StateVector; kwargs...)
DAE(v, u, ϕ, t₀, q₀::State, λ₀::State; kwargs...)
DAE(v, u, ϕ, q₀::State, λ₀::State; kwargs...)

DAE(v, u, ū, ϕ, ψ, t₀, q₀::StateVector, λ₀::StateVector; kwargs...)
DAE(v, u, ū, ϕ, ψ, q₀::StateVector, λ₀::StateVector; kwargs...)
DAE(v, u, ū, ϕ, ψ, t₀, q₀::State, λ₀::State; kwargs...)
DAE(v, u, ū, ϕ, ψ, q₀::State, λ₀::State; kwargs...)
```

### Example

```julia
    function v(t, q, v)
        v[1] = q[1]
        v[2] = q[2]
    end

    function u(t, q, λ, u)
        u[1] = +λ[1]
        u[2] = -λ[1]
    end

    function ϕ(t, q, λ, ϕ)
        ϕ[1] = q[2] - q[1]
    end

    t₀ = 0.
    q₀ = [1., 1.]
    λ₀ = [0.]

    dae = DAE(v, u, ϕ, t₀, q₀, λ₀)

```

"""
struct DAE{dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
           vType <: Function,
           uType <: Function, ūType <: OptionalFunction,
           ϕType <: Function, ψType <: OptionalFunction,
           v̄Type <: Function,
           invType <: OptionalNamedTuple,
           parType <: OptionalNamedTuple,
           perType <: OptionalArray{arrayType}} <: AbstractEquationDAE{dType, tType}

    d::Int
    m::Int
    v::vType
    u::uType
    ū::ūType
    ϕ::ϕType
    ψ::ψType
    v̄::v̄Type
    t₀::tType
    q₀::Vector{arrayType}
    λ₀::Vector{arrayType}
    invariants::invType
    parameters::parType
    periodicity::perType

    function DAE(v::vType, u::uType, ū::ūType, ϕ::ϕType, ψ::ψType, v̄::v̄Type,
                 t₀::tType, q₀::Vector{arrayType}, λ₀::Vector{arrayType},
                 invariants::invType, parameters::parType, periodicity::perType) where {
                dType <: Number, tType <: Number, arrayType <: AbstractArray{dType},
                vType <: Function,
                uType <: Function, ūType <: OptionalFunction,
                ϕType <: Function, ψType <: OptionalFunction,
                v̄Type <: Function,
                invType <: OptionalNamedTuple,
                parType <: OptionalNamedTuple,
                perType <: OptionalArray{arrayType}}

        d = length(q₀[begin])
        m = length(λ₀[begin])

        @assert all(length(q) == d for q in q₀)
        @assert all(length(λ) == m for λ in λ₀)

        new{dType, tType, arrayType, vType, uType, ūType, ϕType, ψType, v̄Type, invType, parType, perType}(d, m, v, u, ū, ϕ, ψ, v̄, t₀, q₀, λ₀, invariants, parameters, periodicity)
    end
end

_DAE(v, u, ū, ϕ, ψ, t₀, q₀, λ₀; v̄=v, invariants=nothing, parameters=nothing, periodicity=nothing) = DAE(v, u, ū, ϕ, ψ, v̄, t₀, q₀, λ₀, invariants, parameters, periodicity)

DAE(v, u, ϕ, t₀, q₀::StateVector, λ₀::StateVector; kwargs...) = _DAE(v, u, nothing, ϕ, nothing, t₀, q₀, λ₀; kwargs...)
DAE(v, u, ϕ, q₀::StateVector, λ₀::StateVector; kwargs...) = DAE(v, u, ϕ, 0.0, q₀, λ₀; kwargs...)
DAE(v, u, ϕ, t₀, q₀::State, λ₀::State; kwargs...) = DAE(v, u, ϕ, t₀, [q₀], [λ₀]; kwargs...)
DAE(v, u, ϕ, q₀::State, λ₀::State; kwargs...) = DAE(v, u, ϕ, 0.0, q₀, λ₀; kwargs...)

DAE(v, u, ū, ϕ, ψ, t₀, q₀::StateVector, λ₀::StateVector; kwargs...) = _DAE(v, u, ū, ϕ, ψ, t₀, q₀, λ₀; kwargs...)
DAE(v, u, ū, ϕ, ψ, q₀::StateVector, λ₀::StateVector; kwargs...) = DAE(v, u, ū, ϕ, ψ, 0.0, q₀, λ₀; kwargs...)
DAE(v, u, ū, ϕ, ψ, t₀, q₀::State, λ₀::State; kwargs...) = DAE(v, u, ū, ϕ, ψ, t₀, [q₀], [λ₀]; kwargs...)
DAE(v, u, ū, ϕ, ψ, q₀::State, λ₀::State; kwargs...) = DAE(v, u, ū, ϕ, ψ, 0.0, q₀, λ₀; kwargs...)

const DAEpsiType{psiT,DT,TT,AT,VT,UT,ŪT,ΦT,V̄T,invT,parT,perT} = DAE{DT,TT,AT,VT,UT,ŪT,ΦT,psiT,V̄T,invT,parT,perT} # type alias for dispatch on secondary constraint type parameter
const DAEinvType{invT,DT,TT,AT,VT,UT,ŪT,ΦT,psiT,V̄T,parT,perT} = DAE{DT,TT,AT,VT,UT,ŪT,ΦT,psiT,V̄T,invT,parT,perT} # type alias for dispatch on invariants type parameter
const DAEparType{parT,DT,TT,AT,VT,UT,ŪT,ΦT,psiT,V̄T,invT,perT} = DAE{DT,TT,AT,VT,UT,ŪT,ΦT,psiT,V̄T,invT,parT,perT} # type alias for dispatch on parameters type parameter
const DAEperType{perT,DT,TT,AT,VT,UT,ŪT,ΦT,psiT,V̄T,invT,parT} = DAE{DT,TT,AT,VT,UT,ŪT,ΦT,psiT,V̄T,invT,parT,perT} # type alias for dispatch on periodicity type parameter

Base.hash(dae::DAE, h::UInt) = hash(dae.d, hash(dae.m, 
        hash(dae.v, hash(dae.u, hash(dae.ū,
        hash(dae.ϕ, hash(dae.ψ, hash(dae.v̄,
        hash(dae.t₀, hash(dae.q₀, hash(dae.λ₀,
        hash(dae.invariants, hash(dae.parameters, hash(dae.periodicity, h))))))))))))))

Base.:(==)(dae1::DAE, dae2::DAE) = (
                                dae1.d == dae2.d
                             && dae1.m == dae2.m
                             && dae1.v == dae2.v
                             && dae1.u == dae2.u
                             && dae1.ū == dae2.ū
                             && dae1.ϕ == dae2.ϕ
                             && dae1.ψ == dae2.ψ
                             && dae1.v̄ == dae2.v̄
                             && dae1.t₀ == dae2.t₀
                             && dae1.q₀ == dae2.q₀
                             && dae1.λ₀ == dae2.λ₀
                             && dae1.invariants  == dae2.invariants
                             && dae1.parameters  == dae2.parameters
                             && dae1.periodicity == dae2.periodicity)

function Base.similar(equ::DAE, t₀::Real, q₀::StateVector, λ₀::StateVector; parameters=equ.parameters)
    @assert all([length(q) == equ.d for q in q₀])
    @assert all([length(λ) == equ.m for λ in λ₀])
    DAE(equ.v, equ.u, equ.ū, equ.ϕ, equ.ψ, t₀, q₀, λ₀; v̄=equ.v̄,
        invariants=equ.invariants, parameters=parameters, periodicity=equ.periodicity)
end

Base.similar(equ::DAE, q₀, λ₀=get_λ₀(q₀, equ.λ₀); kwargs...) = similar(equ, equ.t₀, q₀, λ₀; kwargs...)
Base.similar(equ::DAE, t₀::Real, q₀::State, λ₀::State=get_λ₀(q₀, equ.λ₀); kwargs...) = similar(equ, t₀, [q₀], [λ₀]; kwargs...)

hassecondary(::DAEpsiType{<:Nothing}) = false
hassecondary(::DAEpsiType{<:Function}) = true

hasinvariants(::DAEinvType{<:Nothing}) = false
hasinvariants(::DAEinvType{<:NamedTuple}) = true

hasparameters(::DAEparType{<:Nothing}) = false
hasparameters(::DAEparType{<:NamedTuple}) = true

hasperiodicity(::DAEperType{<:Nothing}) = false
hasperiodicity(::DAEperType{<:AbstractArray}) = true

@inline Base.axes(equation::DAE) = axes(equation.q₀[begin])
@inline Base.ndims(equation::DAE) = equation.d
@inline Common.nsamples(equation::DAE) = length(eachindex(equation.q₀))
@inline Common.nconstraints(equation::DAE) = equation.m

@inline Common.periodicity(equation::DAE) = hasperiodicity(equation) ? equation.periodicity : zero(equation.q₀[begin])
@inline initial_conditions(equation::DAE) = (equation.t₀, equation.q₀, equation.λ₀)

_get_v(equ::DAE) = hasparameters(equ) ? (t,q,v)   -> equ.v(t, q, v, equ.parameters) : equ.v
_get_u(equ::DAE) = hasparameters(equ) ? (t,q,λ,u) -> equ.u(t, q, λ, u, equ.parameters) : equ.u
_get_ū(equ::DAE) = hasparameters(equ) ? (t,q,λ,u) -> equ.ū(t, q, λ, u, equ.parameters) : equ.ū
_get_ϕ(equ::DAE) = hasparameters(equ) ? (t,q,ϕ)   -> equ.ϕ(t, q, ϕ, equ.parameters) : equ.ϕ
_get_ψ(equ::DAE) = hasparameters(equ) ? (t,q,v,ψ) -> equ.ψ(t, q, v, ϕ, equ.parameters) : equ.ψ
_get_v̄(equ::DAE) = hasparameters(equ) ? (t,q,v)   -> equ.v̄(t, q, v, equ.parameters) : equ.v̄

function get_function_tuple(equ::DAE)
    names = (:v,:u,:ϕ,:v̄)
    equs  = (_get_v(equ), _get_u(equ), _get_ϕ(equ), _get_v̄(equ))

    if hassecondary(equ)
        names = (names..., :ū, :ψ)
        equs  = (equs..., _get_ū(equ), _get_ψ(equ))
    end

    NamedTuple{names}(equs)
end
