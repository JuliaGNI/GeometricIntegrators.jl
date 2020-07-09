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

### Fields

* `d`: dimension of dynamical variables ``q`` and ``p`` as well as the vector fields ``v`` and ``f``
* `v`: function computing the vector field ``v``
* `f`: function computing the vector field ``f``
* `h`: function computing the Hamiltonian ``H``
* `t₀`: initial time
* `q₀`: initial condition for dynamical variable ``q``
* `p₀`: initial condition for dynamical variable ``p``

"""
struct HODE{dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
            vType <: Function, fType <: Function, hType <: Function,
            pType <: Union{NamedTuple,Nothing}} <: AbstractEquationPODE{dType, tType}

    d::Int
    v::vType
    f::fType
    h::hType
    t₀::tType
    q₀::Vector{arrayType}
    p₀::Vector{arrayType}
    parameters::pType
    periodicity::Vector{dType}

    function HODE(v::vType, f::fType, h::hType,
                  t₀::tType, q₀::Vector{arrayType}, p₀::Vector{arrayType};
                  parameters::pType=nothing, periodicity=zero(q₀[begin])) where {
                        dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
                        vType <: Function, fType <: Function, hType <: Function,
                        pType <: Union{NamedTuple,Nothing}}

        d = length(q₀[begin])

        @assert length(q₀) == length(p₀)
        @assert all([length(q) == d for q in q₀])
        @assert all([length(p) == d for p in p₀])

        new{dType, tType, arrayType, vType, fType, hType, pType}(d, v, f, h, t₀, q₀, p₀, parameters, periodicity)
    end
end

HODE(v, f, h, q₀::StateVector, p₀::StateVector; kwargs...) = HODE(v, f, h, 0.0, q₀, p₀; kwargs...)
HODE(v, f, h, t₀, q₀::State, p₀::State; kwargs...) = HODE(v, f, h, t₀, [q₀], [p₀]; kwargs...)
HODE(v, f, h, q₀::State, p₀::State; kwargs...) = HODE(v, f, h, 0.0, q₀, p₀; kwargs...)

const HODEPT{PT,DT,TT,AT,VT,FT,HT} = HODE{DT,TT,AT,VT,FT,HT,PT} # type alias for dispatch on parameters type parameter

Base.hash(ode::HODE, h::UInt) = hash(ode.d, hash(ode.v, hash(ode.f, hash(ode.h,
        hash(ode.t₀, hash(ode.q₀, hash(ode.p₀, hash(ode.periodicity, hash(ode.parameters, h)))))))))

Base.:(==)(ode1::HODE, ode2::HODE) = (
                                ode1.d == ode2.d
                             && ode1.v == ode2.v
                             && ode1.f == ode2.f
                             && ode1.h == ode2.h
                             && ode1.t₀ == ode2.t₀
                             && ode1.q₀ == ode2.q₀
                             && ode1.p₀ == ode2.p₀
                             && ode1.parameters == ode2.parameters
                             && ode1.periodicity == ode2.periodicity)

Base.similar(equ::HODE, q₀, p₀; kwargs...) = similar(equ, equ.t₀, q₀, p₀; kwargs...)
Base.similar(equ::HODE, t₀::Real, q₀::State, p₀::State; kwargs...) = similar(equ, t₀, [q₀], [p₀]; kwargs...)

function Base.similar(equ::HODE, t₀::Real, q₀::StateVector, p₀::StateVector;
                      parameters=equ.parameters, periodicity=equ.periodicity)
    @assert all([length(q) == ndims(equ) for q in q₀])
    @assert all([length(p) == ndims(equ) for p in p₀])
    HODE(equ.v, equ.f, equ.h, t₀, q₀, p₀; parameters=parameters, periodicity=periodicity)
end

Base.ndims(equ::HODE) = equ.d
Common.nsamples(equ::HODE) = length(equ.q₀)
Common.periodicity(equ::HODE) = equ.periodicity
initial_conditions(equation::HODE) = (equation.t₀, equation.q₀, equation.p₀)

hashamiltonian(::HODE) = true

hasparameters(::HODEPT{<:Nothing}) = false
hasparameters(::HODEPT{<:NamedTuple}) = true

_get_v(equ::HODE) = hasparameters(equ) ? (t,q,p,v) -> equ.v(t, q, p, v, equ.parameters) : equ.v
_get_f(equ::HODE) = hasparameters(equ) ? (t,q,p,f) -> equ.f(t, q, p, f, equ.parameters) : equ.f
_get_h(equ::HODE) = hasparameters(equ) ? (t,q,p) -> equ.h(t, q, p, equ.parameters) : equ.h

function get_function_tuple(equ::HODE)
    names = (:v,:f,:h)
    equs  = (_get_v(equ), _get_f(equ), _get_h(equ))
    NamedTuple{names}(equs)
end
