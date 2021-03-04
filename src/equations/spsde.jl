@doc raw"""
`SPSDE`: Stratonovich Split Partitioned Stochastic Differential Equation

Defines a partitioned stochastic differential initial value problem
```math
\begin{aligned}
dq (t) &=   v(t, q(t)) \, dt + B(t, q(t)) \circ dW , & q(t_{0}) &= q_{0} , \\
dp (t) &= [ f_1(t, q(t)) + f_2(t, q(t)) ] \, dt + [ G_1(t, q(t)) + G_2(t, q(t)) ] \circ dW , & p(t_{0}) &= p_{0} ,
\end{aligned}
```
with the drift vector fields ``v`` and ``f_i``, diffusion matrices ``B`` and ``G_i``,
initial conditions ``q_{0}`` and ``p_{0}``, the dynamical variables ``(q,p)`` taking
values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}``, and the m-dimensional Wiener process W

### Parameters

* `DT <: Number`: data type
* `TT <: Real`: time step type
* `AT <: AbstractArray{DT}`: array type
* `vType <: Function`: type of `v`
* `f1Type <: Function`: type of `f1`
* `f2Type <: Function`: type of `f2`
* `BType <: Function`: type of `B`
* `G1Type <: Function`: type of `G1`
* `G2Type <: Function`: type of `G2`
* `pType <: Union{NamedTuple,Nothing}`: parameters type

### Fields

* `d`:  dimension of dynamical variable ``q`` and the vector fields ``vi``
* `m`:  dimension of the Wiener process
* `ni`: number of initial conditions
* `ns`: number of sample paths
* `v` :  function computing the drift vector field for the position variable ``q``
* `f1`:  function computing the drift vector field for the momentum variable ``p``
* `f2`:  function computing the drift vector field for the momentum variable ``p``
* `B` :  function computing the d x m diffusion matrix for the position variable ``q``
* `G1`:  function computing the d x m diffusion matrix for the momentum variable ``p``
* `G2`:  function computing the d x m diffusion matrix for the momentum variable ``p``
* `t₀`: initial time
* `q₀`: initial condition for dynamical variable ``q`` (may be a random variable itself)
* `p₀`: initial condition for dynamical variable ``p`` (may be a random variable itself)
* `parameters`: either a `NamedTuple` containing the equations parameters or `nothing`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions

The functions `v`, `f`, `B` and `G`, providing the drift vector fields and diffusion matrices, take four arguments,
`v(t, q, p, v)`, `f(t, q, p, f)`, `B(t, q, p,  B)` and `G(t, q, p, G)`, where `t` is the current time, `(q, p)` is the
current solution vector, and `v`, `f`, `B` and `G` are the variables which hold the result
of evaluating the vector fields ``v``, ``f`` and the matrices ``B``, ``G`` on `t` and `(q,p)`.

### Constructors

```julia
SPSDE(m, ns, v, f1, f2, B, G1, G2, t₀, q₀, p₀; parameters=nothing, periodicity=zero(q₀[begin]))
SPSDE(m, ns, v, f1, f2, B, G1, G2, q₀::StateVector, p₀::StateVector; kwargs...) = SPSDE(m, ns, v, f1, f2, B, G1, G2, 0.0, q₀, p₀; kwargs...)
SPSDE(m, ns, v, f1, f2, B, G1, G2, t₀, q₀::State, p₀::State; kwargs...) = SPSDE(m, ns, v, f1, f2, B, G1, G2, t₀, [q₀], [p₀]; kwargs...)
SPSDE(m, ns, v, f1, f2, B, G1, G2, q₀::State, p₀::State; kwargs...) = SPSDE(m, ns, v, f1, f2, B, G1, G2, 0.0, q₀, p₀; kwargs...)
```

### Example

```julia
    function v(λ, t, q, v)
        v[1] = λ*q[1]
        v[2] = λ*q[2]
    end

    function B(μ, t, q, B)
        B[1] = μ*q[1]
        B[2] = μ*q[2]
    end

    t₀ = 0.
    q₀ = [1., 1.]
    λ  = 2.
    μ  = 1.

    v_sde = (t, q, v) -> v(λ, t, q, v)
    B_sde = (t, q, B) -> B(μ, t, q, B)

    sde = SDE(v_sde, B_sde, t₀, q₀)
```
"""
struct SPSDE{dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
             vType <: Function, f1Type <: Function, f2Type <: Function,
             BType <: Function, G1Type <: Function, G2Type <: Function,
             pType <: Union{NamedTuple,Nothing}} <: AbstractEquationPSDE{dType, tType}
    d::Int
    m::Int
    ns::Int
    v ::vType
    f1::f1Type
    f2::f2Type
    B ::BType
    G1::G1Type
    G2::G2Type
    t₀::tType
    q₀::Vector{arrayType}
    p₀::Vector{arrayType}
    parameters::pType
    periodicity::arrayType

    function SPSDE(m, ns,
                    v::vType, f1::f1Type, f2::f2Type,
                    B::BType, G1::G1Type, G2::G2Type,
                    t₀::tType, q₀::Vector{arrayType}, p₀::Vector{arrayType};
                    parameters::pType=nothing, periodicity=zero(q₀[begin])) where {
                        dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
                        vType <: Function, f1Type <: Function, f2Type <: Function,
                        BType <: Function, G1Type <: Function, G2Type <: Function,
                        pType <: Union{NamedTuple,Nothing}}

        d  = length(q₀[begin])
        ni = length(q₀)

        @assert length(q₀) == length(p₀)
        @assert all(length(q) == d for q in q₀)
        @assert all(length(p) == d for p in p₀)

        @assert ns ≥ 1
        @assert ni == 1 || ns == 1
        # either multiple deterministic initial conditions and one sample path
        # or one deterministic initial condition and multiple sample paths

        new{dType,tType,arrayType,vType,f1Type,f2Type,BType,G1Type,G2Type,pType}(d, m, ns, v, f1, f2, B, G1, G2, t₀, q₀, p₀, parameters, periodicity)
    end
end

SPSDE(m, ns, v, f1, f2, B, G1, G2, q₀::StateVector, p₀::StateVector; kwargs...) = SPSDE(m, ns, v, f1, f2, B, G1, G2, 0.0, q₀, p₀; kwargs...)
SPSDE(m, ns, v, f1, f2, B, G1, G2, t₀, q₀::State, p₀::State; kwargs...) = SPSDE(m, ns, v, f1, f2, B, G1, G2, t₀, [q₀], [p₀]; kwargs...)
SPSDE(m, ns, v, f1, f2, B, G1, G2, q₀::State, p₀::State; kwargs...) = SPSDE(m, ns, v, f1, f2, B, G1, G2, 0.0, q₀, p₀; kwargs...)

const SPSDEPT{PT,DT,TT,AT,VT,F1T,F2T,BT,G1T,G2T} = SPSDE{DT,TT,AT,VT,F1T,F2T,BT,G1T,G2T,PT} # type alias for dispatch on parameters type parameter

Base.hash(sde::SPSDE, h::UInt) = hash(sde.d, hash(sde.m, hash(sde.ns,
                                 hash(sde.v, hash(sde.f1, hash(sde.f2,
                                 hash(sde.B, hash(sde.G1, hash(sde.G2,
                                 hash(sde.t₀, hash(sde.q₀, hash(sde.p₀,
                                 hash(sde.parameters, hash(sde.periodicity, h))))))))))))))

Base.:(==)(sde1::SPSDE, sde2::SPSDE) = (
                                sde1.d == sde2.d
                             && sde1.m == sde2.m
                             && sde1.ns == sde2.ns
                             && sde1.v  == sde2.v
                             && sde1.f1 == sde2.f1
                             && sde1.f2 == sde2.f2
                             && sde1.B  == sde2.B
                             && sde1.G1 == sde2.G1
                             && sde1.G2 == sde2.G2
                             && sde1.t₀ == sde2.t₀
                             && sde1.q₀ == sde2.q₀
                             && sde1.p₀ == sde2.p₀
                             && sde1.parameters  == sde2.parameters
                             && sde1.periodicity == sde2.periodicity)

function Base.similar(equ::SPSDE, t₀::Real, q₀::StateVector, p₀::StateVector, ns::Int=equ.ns)
    @assert all([length(q) == ndims(equ) for q in q₀])
    @assert all([length(p) == ndims(equ) for p in p₀])
    SPSDE(equ.m, ns, equ.v, equ.f1, equ.f2, equ.B, equ.G1, equ.G2, t₀, q₀, p₀; parameters=equ.parameters, periodicity=equ.periodicity)
end

Base.similar(equ::SPSDE, t₀::Real, q₀::State, p₀::State, ns::Int=equ.ns; kwargs...) = similar(equ, t₀, [q₀], [p₀], ns; kwargs...)
Base.similar(equ::SPSDE, q₀::Union{State,StateVector}, p₀::Union{State,StateVector}, ns::Int=equ.ns; kwargs...) = similar(equ, equ.t₀, q₀, p₀, ns; kwargs...)

Base.ndims(equ::SPSDE) = equ.d
Base.axes(equ::SPSDE) = axes(equ.q₀[begin])
Common.nsamples(equ::SPSDE) = length(equ.q₀)
Common.periodicity(equ::SPSDE) = equ.periodicity

initial_conditions(equ::SPSDE) = (equ.t₀, equ.q₀, equ.p₀)

hasparameters(::SPSDEPT{<:Nothing}) = false
hasparameters(::SPSDEPT{<:NamedTuple}) = true

_get_v( equ::SPSDE) = hasparameters(equ) ? (t,q,p,v) -> equ.v( t, q, p, v, equ.parameters) : equ.v
_get_f1(equ::SPSDE) = hasparameters(equ) ? (t,q,p,f) -> equ.f1(t, q, p, f, equ.parameters) : equ.f1
_get_f2(equ::SPSDE) = hasparameters(equ) ? (t,q,p,f) -> equ.f2(t, q, p, f, equ.parameters) : equ.f2
_get_B( equ::SPSDE) = hasparameters(equ) ? (t,q,p,B) -> equ.B( t, q, p, B, equ.parameters) : equ.B
_get_G1(equ::SPSDE) = hasparameters(equ) ? (t,q,p,G) -> equ.G1(t, q, p, G, equ.parameters) : equ.G1
_get_G2(equ::SPSDE) = hasparameters(equ) ? (t,q,p,G) -> equ.G2(t, q, p, G, equ.parameters) : equ.G2


function get_function_tuple(equ::SPSDE)
    names = (:v,:f1,:f2,:B,:G1,:G2)
    equs  = (_get_v(equ), _get_f1(equ), _get_f2(equ), _get_B(equ), _get_G1(equ), _get_G2(equ))

    NamedTuple{names}(equs)
end
