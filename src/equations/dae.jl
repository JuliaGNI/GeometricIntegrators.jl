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

### Fields

* `d`: dimension of dynamical variable ``q`` and the vector field ``v``
* `m`: dimension of algebraic variable ``\lambda`` and the constraint ``\phi``
* `n`: number of initial conditions
* `v`: function computing the vector field
* `u`: function computing the projection
* `ϕ`: algebraic constraint
* `v̄`: function computing an initial guess for the velocity field ``v`` (optional)
* `h`: function computing the Hamiltonian (optional)
* `t₀`: initial time
* `q₀`: initial condition for dynamical variable ``q``
* `λ₀`: initial condition for algebraic variable ``\lambda``

The function `v`, providing the vector field, takes three arguments,
`v(t, q, v)`, the functions `u` and `ϕ`, providing the projection and the
algebraic constraint take four arguments, `u(t, q, λ, u)` and `ϕ(t, q, λ, ϕ)`,
where `t` is the current time, `q` and `λ` are the current solution vectors,
and `v`, `u` and `ϕ` are the vectors which hold the result of evaluating the
vector field ``v``, the projection ``u`` and the algebraic constraint ``\phi``
on `t`, `q` and `λ`.

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
           vType <: Function, uType <: Function,
           ϕType <: Function, v̄Type <: Function, hType <: OptionalFunction,
           pType <: Union{NamedTuple,Nothing}} <: AbstractEquationDAE{dType, tType}

    d::Int
    m::Int
    v::vType
    u::uType
    ϕ::ϕType
    v̄::v̄Type
    h::hType
    t₀::tType
    q₀::Vector{arrayType}
    λ₀::Vector{arrayType}
    parameters::pType
    periodicity::Vector{dType}

    function DAE(v::vType, u::uType, ϕ::ϕType, t₀::tType, q₀::Vector{arrayType}, λ₀::Vector{arrayType};
            v̄::v̄Type=v, h::hType=nothing, parameters::pType=nothing,
            periodicity=zero(q₀[begin])) where {
                dType <: Number, tType <: Number, arrayType <: AbstractArray{dType},
                vType <: Function, uType <: Function, ϕType <: Function,
                v̄Type <: Function, hType <: OptionalFunction, pType <: Union{NamedTuple,Nothing}}

        d = length(q₀[begin])
        m = length(λ₀[begin])

        @assert all(length(q) == d for q in q₀)
        @assert all(length(λ) == m for λ in λ₀)

        new{dType, tType, arrayType, vType, uType, ϕType, hType, pType}(d, m, v, u, ϕ, v̄, h, t₀, q₀, λ₀, parameters, periodicity)
    end
end

DAE(v, u, ϕ, q₀::StateVector, λ₀::StateVector; kwargs...) = DAE(v, u, ϕ, 0.0, q₀, λ₀; kwargs...)
DAE(v, u, ϕ, t₀, q₀::State, λ₀::State; kwargs...) = DAE(v, u, ϕ, t₀, [q₀], [λ₀]; kwargs...)
DAE(v, u, ϕ, q₀::State, λ₀::State; kwargs...) = DAE(v, u, ϕ, 0.0, q₀, λ₀; kwargs...)

const DAEHT{HT,DT,TT,AT,VT,UT,ΦT,PT} = DAE{DT,TT,AT,VT,UT,ΦT,HT,PT} # type alias for dispatch on Hamiltonian type parameter
const DAEPT{PT,DT,TT,AT,VT,UT,ΦT,HT} = DAE{DT,TT,AT,VT,UT,ΦT,HT,PT} # type alias for dispatch on parameters type parameter

Base.hash(dae::DAE, h::UInt) = hash(dae.d, hash(dae.m, hash(dae.v,
        hash(dae.u, hash(dae.ϕ, hash(dae.v̄, hash(dae.h, hash(dae.t₀, hash(dae.q₀, hash(dae.λ₀,
        hash(dae.periodicity, hash(dae.parameters, h))))))))))))

Base.:(==)(dae1::DAE, dae2::DAE) = (
                                dae1.d == dae2.d
                             && dae1.m == dae2.m
                             && dae1.v == dae2.v
                             && dae1.u == dae2.u
                             && dae1.ϕ == dae2.ϕ
                             && dae1.v̄ == dae2.v̄
                             && dae1.h == dae2.h
                             && dae1.t₀ == dae2.t₀
                             && dae1.q₀ == dae2.q₀
                             && dae1.λ₀ == dae2.λ₀
                             && dae1.parameters == dae1.parameters
                             && dae1.periodicity == dae1.periodicity)

Base.similar(equ::DAE, q₀, λ₀=get_λ₀(q₀, equ.λ₀); kwargs...) = similar(equ, equ.t₀, q₀, λ₀; kwargs...)
Base.similar(equ::DAE, t₀::Real, q₀::State, λ₀::State=get_λ₀(q₀, equ.λ₀); kwargs...) = similar(equ, t₀, [q₀], [λ₀]; kwargs...)

function Base.similar(equ::DAE, t₀::Real, q₀::StateVector, λ₀::StateVector;
                      v̄=equ.v̄, h=equ.h, parameters=equ.parameters, periodicity=equ.periodicity)
    @assert all([length(q) == equ.d for q in q₀])
    @assert all([length(λ) == equ.m for λ in λ₀])
    DAE(equ.v, equ.u, equ.ϕ, t₀, q₀, λ₀; v̄=v̄, h=h, parameters=parameters, periodicity=periodicity)
end

@inline Base.ndims(equation::DAE) = equation.d
@inline Common.nsamples(equation::DAE) = length(eachindex(equation.q₀))
@inline Common.nconstraints(equation::DAE) = equation.m
@inline Common.periodicity(equation::DAE) = equation.periodicity

initial_conditions(equation::DAE) = (equation.t₀, equation.q₀, equation.λ₀)

hashamiltonian(::DAEHT{<:Nothing}) = false
hashamiltonian(::DAEHT{<:Function}) = true

hasparameters(::DAEPT{<:Nothing}) = false
hasparameters(::DAEPT{<:NamedTuple}) = true

_get_v(equ::DAE) = hasparameters(equ) ? (t,q,v) -> equ.v(t, q, v, equ.parameters) : equ.v
_get_u(equ::DAE) = hasparameters(equ) ? (t,q,λ,u) -> equ.u(t, q, λ, u, equ.parameters) : equ.u
_get_ϕ(equ::DAE) = hasparameters(equ) ? (t,q,ϕ) -> equ.ϕ(t, q, ϕ, equ.parameters) : equ.ϕ
_get_h(equ::DAE) = hasparameters(equ) ? (t,q) -> equ.h(t, q, equ.parameters) : equ.h
_get_v̄(equ::DAE) = hasparameters(equ) ? (t,q,v) -> equ.v̄(t, q, v, equ.parameters) : equ.v̄

function get_function_tuple(equ::DAE)
    names = (:v,:u,:ϕ,:v̄)
    equs  = (_get_v(equ), _get_u(equ), _get_ϕ(equ), _get_v̄(equ))

    if hashamiltonian(equ)
        names = (names..., :h)
        equs  = (equs..., _get_h(equ))
    end

    NamedTuple{names}(equs)
end
