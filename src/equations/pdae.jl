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

### Fields

* `d`: dimension of dynamical variables ``q`` and ``p`` as well as the vector fields ``f`` and ``p``
* `m`: dimension of algebraic variable ``\lambda`` and the constraint ``\phi``
* `v`: function computing the vector field ``v``
* `f`: function computing the vector field ``f``
* `u`: function computing the projection for ``q``
* `g`: function computing the projection for ``p``
* `ϕ`: algebraic constraint
* `v̄`: function computing an initial guess for the velocity field ``v`` (optional)
* `f̄`: function computing an initial guess for the force field ``f`` (optional)
* `h`: function computing the Hamiltonian (optional)
* `t₀`: initial time
* `q₀`: initial condition for dynamical variable ``q``
* `p₀`: initial condition for dynamical variable ``p``
* `λ₀`: initial condition for algebraic variable ``\lambda``

"""
struct PDAE{dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
            vType <: Function, fType <: Function,
            uType <: Function, gType <: Function,
            ϕType <: Function,
            v̄Type <: Function, f̄Type <: Function,
            hType <: OptionalFunction,
            pType <: Union{NamedTuple,Nothing}} <: AbstractEquationPDAE{dType, tType}

    d::Int
    m::Int
    v::vType
    f::fType
    u::uType
    g::gType
    ϕ::ϕType
    v̄::v̄Type
    f̄::f̄Type
    h::hType
    t₀::tType
    q₀::Vector{arrayType}
    p₀::Vector{arrayType}
    λ₀::Vector{arrayType}
    parameters::pType
    periodicity::arrayType

    function PDAE(v::vType, f::fType, u::uType, g::gType, ϕ::ϕType,
                  t₀::tType, q₀::Vector{arrayType}, p₀::Vector{arrayType}, λ₀::Vector{arrayType};
                  v̄::v̄Type=v, f̄::f̄Type=f, h::hType=nothing, parameters::pType=nothing,
                  periodicity=zero(q₀[begin])) where {
                        dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
                        vType <: Function, fType <: Function,
                        uType <: Function, gType <: Function,
                        ϕType <: Function,
                        v̄Type <: Function, f̄Type <: Function,
                        hType <: OptionalFunction,
                        pType <: Union{NamedTuple,Nothing}}

        d = length(q₀[begin])
        m = length(λ₀[begin])

        @assert 2d ≥ m

        @assert length(q₀) == length(p₀) == length(λ₀)

        @assert all(length(q) == d for q in q₀)
        @assert all(length(p) == d for p in p₀)
        @assert all(length(λ) == m for λ in λ₀)

        @assert all([ndims(q) == ndims(p) == ndims(λ) for (q,p,λ) in zip(q₀,p₀,λ₀)])

        new{dType, tType, arrayType, vType, fType, uType, gType, v̄Type, f̄Type, hType, pType}(d, m, v, f, u, g, ϕ, v̄, f̄, h, t₀, q₀, p₀, λ₀, parameters, periodicity)
    end
end

PDAE(v, f, u, g, ϕ, q₀::StateVector, p₀::StateVector, λ₀::StateVector; kwargs...) = PDAE(v, f, u, g, ϕ, 0.0, q₀, p₀, λ₀; kwargs...)
PDAE(v, f, u, g, ϕ, t₀, q₀::State, p₀::State, λ₀::State; kwargs...) = PDAE(v, f, u, g, ϕ, t₀, [q₀], [p₀], [λ₀]; kwargs...)
PDAE(v, f, u, g, ϕ, q₀::State, p₀::State, λ₀::State; kwargs...) = PDAE(v, f, u, g, ϕ, 0.0, q₀, p₀, λ₀; kwargs...)

const PDAEHT{HT,DT,TT,AT,VT,FT,UT,GT,ΦT,PT} = PDAE{DT,TT,AT,VT,FT,UT,GT,ΦT,HT,PT} # type alias for dispatch on Hamiltonian type parameter
const PDAEPT{PT,DT,TT,AT,VT,FT,UT,GT,ΦT,HT} = PDAE{DT,TT,AT,VT,FT,UT,GT,ΦT,HT,PT} # type alias for dispatch on parameters type parameter

Base.hash(dae::PDAE, h::UInt) = hash(dae.d, hash(dae.m,
        hash(dae.v, hash(dae.f, hash(dae.u, hash(dae.g,
        hash(dae.v̄, hash(dae.f̄, hash(dae.h,
        hash(dae.t₀, hash(dae.q₀, hash(dae.p₀, hash(dae.λ₀,
        hash(dae.periodicity, hash(dae.parameters, h)))))))))))))))

Base.:(==)(dae1::PDAE, dae2::PDAE) = (
                                dae1.d == dae2.d
                             && dae1.m == dae2.m
                             && dae1.v == dae2.v
                             && dae1.f == dae2.f
                             && dae1.u == dae2.u
                             && dae1.g == dae2.g
                             && dae1.ϕ == dae2.ϕ
                             && dae1.v̄ == dae2.v̄
                             && dae1.f̄ == dae2.f̄
                             && dae1.h == dae2.h
                             && dae1.t₀ == dae2.t₀
                             && dae1.q₀ == dae2.q₀
                             && dae1.p₀ == dae2.p₀
                             && dae1.λ₀ == dae2.λ₀
                             && dae1.parameters == dae1.parameters
                             && dae1.periodicity == dae1.periodicity)

Base.similar(equ::PDAE, q₀, p₀, λ₀=get_λ₀(q₀, equ.λ₀); kwargs...) = similar(equ, equ.t₀, q₀, p₀, λ₀; kwargs...)
Base.similar(equ::PDAE, t₀::Real, q₀::State, p₀::State, λ₀::State=get_λ₀(q₀, equ.λ₀); kwargs...) = similar(equ, t₀, [q₀], [p₀], [λ₀]; kwargs...)

function Base.similar(equ::PDAE, t₀::Real, q₀::StateVector, p₀::StateVector, λ₀::StateVector;
                      v̄=equ.v̄, f̄=equ.f̄, h=equ.h, parameters=equ.parameters, periodicity=equ.periodicity)
    @assert all([length(q) == equ.d for q in q₀])
    @assert all([length(p) == equ.d for p in p₀])
    @assert all([length(λ) == equ.m for λ in λ₀])
    PDAE(equ.v, equ.f, equ.u, equ.g, equ.ϕ, t₀, q₀, p₀, λ₀; v̄=v̄, f̄=f̄, h=h, parameters=parameters, periodicity=periodicity)
end

@inline Base.ndims(equation::PDAE) = equation.d
@inline Base.axes(equation::PDAE) = axes(equation.q₀[begin])
@inline Common.nsamples(equation::PDAE) = length(eachindex(equation.q₀))
@inline Common.nconstraints(equation::PDAE) = equation.m
@inline Common.periodicity(equation::PDAE) = equation.periodicity

initial_conditions(equation::PDAE) = (equation.t₀, equation.q₀, equation.p₀, equation.λ₀)

hashamiltonian(::PDAEHT{<:Nothing}) = false
hashamiltonian(::PDAEHT{<:Function}) = true

hasparameters(::PDAEPT{<:Nothing}) = false
hasparameters(::PDAEPT{<:NamedTuple}) = true

_get_v(equ::PDAE) = hasparameters(equ) ? (t,q,p,v) -> equ.v(t, q, p, v, equ.parameters) : equ.v
_get_f(equ::PDAE) = hasparameters(equ) ? (t,q,p,f) -> equ.f(t, q, p, f, equ.parameters) : equ.f
_get_u(equ::PDAE) = hasparameters(equ) ? (t,q,p,λ,u) -> equ.u(t, q, p, λ, u, equ.parameters) : equ.u
_get_g(equ::PDAE) = hasparameters(equ) ? (t,q,p,λ,g) -> equ.g(t, q, p, λ, g, equ.parameters) : equ.g
_get_ϕ(equ::PDAE) = hasparameters(equ) ? (t,q,p,ϕ) -> equ.ϕ(t, q, p, ϕ, equ.parameters) : equ.ϕ
_get_v̄(equ::PDAE) = hasparameters(equ) ? (t,q,p,v) -> equ.v̄(t, q, p, v, equ.parameters) : equ.v̄
_get_f̄(equ::PDAE) = hasparameters(equ) ? (t,q,p,f) -> equ.f̄(t, q, p, f, equ.parameters) : equ.f̄
_get_h(equ::PDAE) = hasparameters(equ) ? (t,q,p) -> equ.h(t, q, p, equ.parameters) : equ.h


function get_function_tuple(equ::PDAE)
    names = (:v,:f,:u,:g,:ϕ,:v̄,:f̄)
    equs  = (_get_v(equ), _get_f(equ), _get_u(equ), _get_g(equ), _get_ϕ(equ), _get_v̄(equ), _get_f̄(equ))

    if hashamiltonian(equ)
        names = (names..., :h)
        equs  = (equs..., _get_h(equ))
    end

    NamedTuple{names}(equs)
end
