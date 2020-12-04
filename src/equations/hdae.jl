@doc raw"""
`HDAE`: Hamiltonian Differential Algebraic Equation *EXPERIMENTAL*

Defines a Hamiltonian differential algebraic initial value problem, that is
a canonical Hamiltonian system of equations subject to Dirac constraints,
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t), p(t)) + u(t, q(t), p(t), \lambda(t)) + \bar{g}(t, q(t), p(t), \lambda(t), \gamma(t)) , & q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), p(t)) + g(t, q(t), p(t), \lambda(t)) + \bar{f}(t, q(t), p(t), \lambda(t), \gamma(t)) , & p(t_{0}) &= p_{0} , \\
0 &= \phi (t, q(t), p(t)) , \\
0 &= \psi (t, q(t), p(t), \dot{q}(t), \dot{p}(t)) ,
\end{aligned}
```
with vector fields ``v``, ``u``, ``\bar{u}`` and ``f``, ``g``, ``\bar{g}``,
primary constraint ``\phi(q,p)=0`` and secondary constraint ``\psi(q,p,\lambda)=0``,
initial conditions ``(q_{0}, p_{0})``, the dynamical variables ``(q,p)``
taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}`` and
the algebraic variables ``(\lambda, \gamma)`` taking values in
``\mathbb{R}^{n} \times \mathbb{R}^{d}``.

### Fields

* `d`: dimension of dynamical variables ``q`` and ``p`` as well as the vector fields ``v`` and ``f``
* `m`: dimension of algebraic variables ``\lambda`` and ``\gamma`` and the constraints ``\phi`` and ``\psi``
* `v`: function computing the Hamiltonian vector field ``v``
* `f`: function computing the Hamiltonian vector field ``f``
* `u`: function computing the primary projection field ``u``
* `g`: function computing the primary projection field ``g``
* `u̅`: function computing the secondary projection field ``\bar{u}``
* `g̅`: function computing the secondary projection field ``\bar{g}``
* `ϕ`: primary constraints
* `ψ`: secondary constraints
* `v̄`: function computing an initial guess for the velocity field ``v``` (optional)
* `f̄`: function computing an initial guess for the force field ``f`` (optional)
* `h`: function computing the Hamiltonian ``H``
* `t₀`: initial time
* `q₀`: initial condition for dynamical variable ``q``
* `p₀`: initial condition for dynamical variable ``p``
* `λ₀`: initial condition for algebraic variable ``λ``
"""
struct HDAE{dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
            vType <: Function, fType <: Function,
            uType <: Function, gType <: Function,
            u̅Type <: Function, g̅Type <: Function,
            ϕType <: Function, ψType <: Function,
            hType <: Function,
            v̄Type <: Function, f̄Type <: Function,
            pType <: Union{NamedTuple,Nothing}} <: AbstractEquationPDAE{dType, tType}
    d::Int
    m::Int
    v::vType
    f::fType
    u::uType
    g::gType
    u̅::u̅Type
    g̅::g̅Type
    ϕ::ϕType
    ψ::ψType
    h::hType
    v̄::v̄Type
    f̄::f̄Type
    t₀::tType
    q₀::Vector{arrayType}
    p₀::Vector{arrayType}
    λ₀::Vector{arrayType}
    parameters::pType
    periodicity::arrayType

    function HDAE(v::vType, f::fType, u::uType, g::gType,
                  u̅::u̅Type, g̅::g̅Type, ϕ::ϕType, ψ::ψType, h::hType,
                  t₀::tType, q₀::Vector{arrayType}, p₀::Vector{arrayType}, λ₀::Vector{arrayType};
                  v̄::v̄Type=v, f̄::f̄Type=f, parameters::pType=nothing, periodicity=zero(q₀[begin])) where {
                        dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
                        vType <: Function, fType <: Function,
                        uType <: Function, gType <: Function,
                        u̅Type <: Function, g̅Type <: Function,
                        ϕType <: Function, ψType <: Function,
                        hType <: Function,
                        v̄Type <: Function, f̄Type <: Function,
                        pType <: Union{NamedTuple,Nothing}}

        d = length(q₀[begin])
        m = length(λ₀[begin])

        @assert 2d ≥ m

        @assert length(q₀) == length(p₀) == length(λ₀)

        @assert all(length(q) == d for q in q₀)
        @assert all(length(p) == d for p in p₀)
        @assert all(length(λ) == m for λ in λ₀)

        @assert all([ndims(q) == ndims(p) == ndims(λ) for (q,p,λ) in zip(q₀,p₀,λ₀)])

        new{dType, tType, arrayType, vType, fType, uType, gType, u̅Type, g̅Type, ϕType, ψType, hType, v̄Type, f̄Type, pType}(
                d, m, v, f, u, g, u̅, g̅, ϕ, ψ, h, v̄, f̄, t₀, q₀, p₀, λ₀, parameters, periodicity)
    end
end

HDAE(v, f, u, g, u̅, g̅, ϕ, ψ, h, q₀::StateVector, p₀::StateVector, λ₀::StateVector; kwargs...) = HDAE(v, f, u, g, u̅, g̅, ϕ, ψ, h, 0.0, q₀, p₀, λ₀; kwargs...)
HDAE(v, f, u, g, u̅, g̅, ϕ, ψ, h, t₀, q₀::State, p₀::State, λ₀::State; kwargs...) = HDAE(v, f, u, g, u̅, g̅, ϕ, ψ, h, t₀, [q₀], [p₀], [λ₀]; kwargs...)
HDAE(v, f, u, g, u̅, g̅, ϕ, ψ, h, q₀::State, p₀::State, λ₀::State; kwargs...) = HDAE(v, f, u, g, u̅, g̅, ϕ, ψ, h, 0.0, q₀, p₀, λ₀; kwargs...)

const HDAEPT{PT,DT,TT,AT,VT,FT,UT,GT,U̅T,G̅T,ΦT,ΨT,HT,V̄T,F̄T} = HDAE{DT,TT,AT,VT,FT,UT,GT,U̅T,G̅T,ΦT,ΨT,HT,V̄T,F̄T,PT} # type alias for dispatch on parameters type parameter

Base.hash(dae::HDAE, h::UInt) = hash(dae.d, hash(dae.m,
                        hash(dae.v, hash(dae.f, hash(dae.u, hash(dae.g,
                        hash(dae.u̅, hash(dae.g̅, hash(dae.ϕ, hash(dae.ψ,
                        hash(dae.h, hash(dae.v̄, hash(dae.f̄,
                        hash(dae.t₀, hash(dae.q₀, hash(dae.p₀, h))))))))))))))))

Base.:(==)(dae1::HDAE, dae2::HDAE) = (
                                dae1.d == dae2.d
                             && dae1.m == dae2.m
                             && dae1.v == dae2.v
                             && dae1.f == dae2.f
                             && dae1.u == dae2.u
                             && dae1.g == dae2.g
                             && dae1.u̅ == dae2.u̅
                             && dae1.g̅ == dae2.g̅
                             && dae1.ϕ == dae2.ϕ
                             && dae1.ψ == dae2.ψ
                             && dae1.h == dae2.h
                             && dae1.v̄ == dae2.v̄
                             && dae1.f̄ == dae2.f̄
                             && dae1.t₀ == dae2.t₀
                             && dae1.q₀ == dae2.q₀
                             && dae1.p₀ == dae2.p₀)

Base.similar(equ::HDAE, q₀, p₀, λ₀=get_λ₀(q₀, equ.λ₀); kwargs...) = similar(equ, equ.t₀, q₀, p₀, λ₀; kwargs...)
Base.similar(equ::HDAE, t₀::Real, q₀::State, p₀::State, λ₀::State=get_λ₀(q₀, equ.λ₀); kwargs...) = similar(equ, t₀, [q₀], [p₀], [λ₀]; kwargs...)

function Base.similar(equ::HDAE, t₀::Real, q₀::StateVector, p₀::StateVector, λ₀::StateVector;
                      v̄=equ.v̄, f̄=equ.f̄, parameters=equ.parameters, periodicity=equ.periodicity)
    @assert all([length(q) == equ.d for q in q₀])
    @assert all([length(p) == equ.d for p in p₀])
    @assert all([length(λ) == equ.m for λ in λ₀])
    HDAE(equ.v, equ.f, equ.u, equ.g, equ.u̅, equ.g̅, equ.ϕ, equ.ψ, equ.h, t₀, q₀, p₀, λ₀;
         v̄=v̄, f̄=f̄, parameters=parameters, periodicity=periodicity)
end


@inline Base.ndims(equation::HDAE) = equation.d
@inline Base.axes(equation::HDAE) = axes(equation.q₀[begin])
@inline Common.nsamples(equation::HDAE) = length(eachindex(equation.q₀))
@inline Common.nconstraints(equation::HDAE) = equation.m
@inline Common.periodicity(equation::HDAE) = equation.periodicity

initial_conditions(equation::HDAE) = (equation.t₀, equation.q₀, equation.p₀, equation.λ₀)

hashamiltonian(::HDAE) = true

hasparameters(::HDAEPT{<:Nothing}) = false
hasparameters(::HDAEPT{<:NamedTuple}) = true

_get_v(equ::HDAE) = hasparameters(equ) ? (t,q,p,v) -> equ.v(t, q, p, v, equ.parameters) : equ.v
_get_f(equ::HDAE) = hasparameters(equ) ? (t,q,p,f) -> equ.f(t, q, p, f, equ.parameters) : equ.f
_get_u(equ::HDAE) = hasparameters(equ) ? (t,q,p,λ,u) -> equ.u(t, q, p, λ, u, equ.parameters) : equ.u
_get_g(equ::HDAE) = hasparameters(equ) ? (t,q,p,λ,g) -> equ.g(t, q, p, λ, g, equ.parameters) : equ.g
_get_u̅(equ::HDAE) = hasparameters(equ) ? (t,q,p,λ,u̅) -> equ.u̅(t, q, p, λ, u̅, equ.parameters) : equ.u̅
_get_g̅(equ::HDAE) = hasparameters(equ) ? (t,q,p,λ,g̅) -> equ.g̅(t, q, p, λ, g̅, equ.parameters) : equ.g̅
_get_ϕ(equ::HDAE) = hasparameters(equ) ? (t,q,p,ϕ) -> equ.ϕ(t, q, p, ϕ, equ.parameters) : equ.ϕ
_get_ψ(equ::HDAE) = hasparameters(equ) ? (t,q,p,v,f,ψ) -> equ.ψ(t, q, p, v, f, ψ, equ.parameters) : equ.ψ
_get_h(equ::HDAE) = hasparameters(equ) ? (t,q,p) -> equ.h(t, q, p, equ.parameters) : equ.h
_get_v̄(equ::HDAE) = hasparameters(equ) ? (t,q,p,v) -> equ.v̄(t, q, p, v, equ.parameters) : equ.v̄
_get_f̄(equ::HDAE) = hasparameters(equ) ? (t,q,p,f) -> equ.f̄(t, q, p, f, equ.parameters) : equ.f̄

function get_function_tuple(equ::HDAE)
    NamedTuple{(:v, :f, :u, :g, :u̅, :g̅, :ϕ, :ψ, :h, :v̄, :f̄)}((
        _get_v(equ), _get_f(equ), _get_u(equ), _get_g(equ),
        _get_u̅(equ), _get_g̅(equ), _get_ϕ(equ), _get_ψ(equ),
        _get_h(equ), _get_v̄(equ), _get_f̄(equ)))
end
