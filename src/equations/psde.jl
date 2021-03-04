@doc raw"""
`PSDE`: Stratonovich Partitioned Stochastic Differential Equation

Defines a partitioned stochastic differential initial value problem
```math
\begin{aligned}
dq (t) &= v(t, q(t)) \, dt + B(t, q(t)) \circ dW , & q(t_{0}) &= q_{0} , \\
dp (t) &= f(t, q(t)) \, dt + G(t, q(t)) \circ dW , & p(t_{0}) &= p_{0}
\end{aligned}
```
with the drift vector fields ``v`` and ``f``, diffusion matrices ``B`` and ``G``,
initial conditions ``q_{0}`` and ``p_{0}``, the dynamical variables ``(q,p)`` taking
values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}``, and the m-dimensional Wiener process W

### Parameters

* `DT <: Number`: data type
* `TT <: Real`: time step type
* `AT <: AbstractArray{DT}`: array type
* `vType <: Function`: type of `v`
* `fType <: Function`: type of `f`
* `BType <: Function`: type of `B`
* `GType <: Function`: type of `G`
* `pType <: Union{NamedTuple,Nothing}`: parameters type

### Fields

* `d`:  dimension of dynamical variable ``q`` and the vector field ``v``
* `m`:  dimension of the Wiener process
* `ns`: number of sample paths
* `v`:  function computing the drift vector field for the position variable ``q``
* `f`:  function computing the drift vector field for the momentum variable ``p``
* `B`:  function computing the d x m diffusion matrix for the position variable ``q``
* `G`:  function computing the d x m diffusion matrix for the momentum variable ``p``
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
PSDE(m, ns, v, f, B, G, t₀, q₀, p₀; parameters=nothing, periodicity=zero(q₀[begin]))
PSDE(m, ns, v, f, B, G, q₀::StateVector, p₀::StateVector; kwargs...) = PSDE(m, ns, v, f, B, G, 0.0, q₀, p₀; kwargs...)
PSDE(m, ns, v, f, B, G, t₀, q₀::State, p₀::State; kwargs...) = PSDE(m, ns, v, f, B, G, t₀, [q₀], [p₀]; kwargs...)
PSDE(m, ns, v, f, B, G, q₀::State, p₀::State; kwargs...) = PSDE(m, ns, v, f, B, G, 0.0, q₀, p₀; kwargs...)
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
struct PSDE{dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
            vType <: Function, fType <: Function,
            BType <: Function, GType <: Function,
            pType <: Union{NamedTuple,Nothing}} <: AbstractEquationPSDE{dType, tType}

    d::Int
    m::Int
    ns::Int
    v::vType
    f::fType
    B::BType
    G::GType
    t₀::tType
    q₀::Vector{arrayType}
    p₀::Vector{arrayType}
    parameters::pType
    periodicity::arrayType

    function PSDE(m, ns, v::vType, f::fType, B::BType, G::GType,
                  t₀::tType, q₀::Vector{arrayType}, p₀::Vector{arrayType};
                  parameters::pType=nothing, periodicity=zero(q₀[begin])) where {
                        dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
                        vType <: Function, fType <: Function,
                        BType <: Function, GType <: Function,
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

        new{dType, tType, arrayType, vType, fType, BType, GType, pType}(d, m, ns, v, f, B, G, t₀, q₀, p₀, parameters, periodicity)
    end
end

PSDE(m, ns, v, f, B, G, q₀::StateVector, p₀::StateVector; kwargs...) = PSDE(m, ns, v, f, B, G, 0.0, q₀, p₀; kwargs...)
PSDE(m, ns, v, f, B, G, t₀, q₀::State, p₀::State; kwargs...) = PSDE(m, ns, v, f, B, G, t₀, [q₀], [p₀]; kwargs...)
PSDE(m, ns, v, f, B, G, q₀::State, p₀::State; kwargs...) = PSDE(m, ns, v, f, B, G, 0.0, q₀, p₀; kwargs...)

# const PSDEHT{HT,DT,TT,AT,VT,FT,BT,GT,PT} = PSDE{DT,TT,AT,VT,FT,BT,GT,HT,PT} # type alias for dispatch on Hamiltonian type parameter
# const PSDEPT{PT,DT,TT,AT,VT,FT,BT,GT,HT} = PSDE{DT,TT,AT,VT,FT,BT,GT,HT,PT} # type alias for dispatch on parameters type parameter
const PSDEPT{PT,DT,TT,AT,VT,FT,BT,GT} = PSDE{DT,TT,AT,VT,FT,BT,GT,PT} # type alias for dispatch on parameters type parameter


Base.hash(sde::PSDE, h::UInt) = hash(sde.d, hash(sde.m, hash(sde.ns,
                                hash(sde.v, hash(sde.f, hash(sde.B, hash(sde.G,
                                hash(sde.t₀, hash(sde.q₀, hash(sde.p₀,
                                hash(sde.parameters, hash(sde.periodicity, h))))))))))))

Base.:(==)(sde1::PSDE, sde2::PSDE) = (
                                sde1.d == sde2.d
                             && sde1.m == sde2.m
                             && sde1.ns == sde2.ns
                             && sde1.v == sde2.v
                             && sde1.f == sde2.f
                             && sde1.B == sde2.B
                             && sde1.G == sde2.G
                             && sde1.t₀ == sde2.t₀
                             && sde1.q₀ == sde2.q₀
                             && sde1.p₀ == sde2.p₀
                             && sde1.parameters  == sde2.parameters
                             && sde1.periodicity == sde2.periodicity)

function Base.similar(equ::PSDE, t₀::Real, q₀::StateVector, p₀::StateVector, ns::Int=equ.ns)
    @assert all([length(q) == ndims(equ) for q in q₀])
    @assert all([length(p) == ndims(equ) for p in p₀])
    PSDE(equ.m, ns, equ.v, equ.f, equ.B, equ.G, t₀, q₀, p₀; parameters=equ.parameters, periodicity=equ.periodicity)
end

Base.similar(equ::PSDE, t₀::Real, q₀::State, p₀::State, ns::Int=equ.ns; kwargs...) = similar(equ, t₀, [q₀], [p₀], ns; kwargs...)
Base.similar(equ::PSDE, q₀::Union{State,StateVector}, p₀::Union{State,StateVector}, ns::Int=equ.ns; kwargs...) = similar(equ, equ.t₀, q₀, p₀, ns; kwargs...)

Base.ndims(equ::PSDE) = equ.d
Base.axes(equ::PSDE) = axes(equ.q₀[begin])
Common.nsamples(equ::PSDE) = length(equ.q₀)
Common.periodicity(equ::PSDE) = equ.periodicity

initial_conditions(equ::PSDE) = (equ.t₀, equ.q₀, equ.p₀)

# hashamiltonian(::PSDEHT{<:Nothing}) = false
# hashamiltonian(::PSDEHT{<:Function}) = true

hasparameters(::PSDEPT{<:Nothing}) = false
hasparameters(::PSDEPT{<:NamedTuple}) = true

_get_v(equ::PSDE) = hasparameters(equ) ? (t,q,p,v) -> equ.v(t, q, p, v, equ.parameters) : equ.v
_get_f(equ::PSDE) = hasparameters(equ) ? (t,q,p,f) -> equ.f(t, q, p, f, equ.parameters) : equ.f
_get_B(equ::PSDE) = hasparameters(equ) ? (t,q,p,B) -> equ.B(t, q, p, B, equ.parameters) : equ.B
_get_G(equ::PSDE) = hasparameters(equ) ? (t,q,p,G) -> equ.G(t, q, p, G, equ.parameters) : equ.G
# _get_h(equ::PSDE) = hasparameters(equ) ? (t,q,p) -> equ.h(t, q, p, equ.parameters) : equ.h


function get_function_tuple(equ::PSDE)
    names = (:v,:f,:B,:G)
    equs  = (_get_v(equ), _get_f(equ), _get_B(equ), _get_G(equ))

    # if hashamiltonian(equ)
    #     names = (names..., :h)
    #     equs  = (equs..., _get_h(equ))
    # end

    NamedTuple{names}(equs)
end
