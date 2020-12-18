@doc raw"""
`SDE`: Stratonovich Stochastic Differential Equation

Defines a stochastic differential initial value problem
```math
\begin{aligned}
dq (t) &= v(t, q(t)) \, dt + B(t, q(t)) \circ dW , & q(t_{0}) &= q_{0} ,
\end{aligned}
```
with drift vector field ``v``, diffusion matrix ``B``,
initial conditions ``q_{0}``, the dynamical variable ``q``
taking values in ``\mathbb{R}^{d}``, and the m-dimensional Wiener process W

### Fields

* `d`:  dimension of dynamical variable ``q`` and the vector field ``v``
* `m`:  dimension of the Wiener process
* `ns`: number of sample paths
* `v`:  function computing the deterministic vector field
* `B`:  function computing the d x m diffusion matrix
* `t₀`: initial time
* `q₀`: initial condition for dynamical variable ``q`` (may be a random variable itself)

### Parameters

* `N`: dimension of nitial condition array: N=1 - single, N=2 - multiple


The functions `v` and `B`, providing the drift vector field and diffusion matrix,
`v(t, q, v)` and `B(t, q, B, col=0)`, where `t` is the current time, `q` is the
current solution vector, and `v` and `B` are the variables which hold the result
of evaluating the vector field ``v`` and the matrix ``B`` on `t` and `q` (if col==0),
or the column col of the matrix B (if col>0).

### Example

```julia
    function v(t, q, v, p)
        λ = p[:λ]
        v[1] = λ*q[1]
        v[2] = λ*q[2]
    end

    function B(t, q, B, p, col=0)
        μ = p[:μ]
        if col==0 #whole matrix
            B[1,1] = μ*q[1]
            B[2,1] = μ*q[2]
        elseif col==1
            #just first column
        end
    end

    t₀ = 0.
    q₀ = [1., 1.]
    λ  = 2.
    μ  = 1.
    p = (λ=λ, μ=μ)

    sde = SDE(v, B, t₀, q₀; parameters=p)
```
"""
struct SDE{dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
           vType <: Function, BType <: Function,
           pType <: Union{NamedTuple,Nothing}} <: AbstractEquationSDE{dType, tType}

    d::Int
    m::Int
    ns::Int
    v::vType
    B::BType
    t₀::tType
    q₀::Vector{arrayType}
    parameters::pType
    periodicity::arrayType

    function SDE(m, ns, v::vType, B::BType, t₀::tType, q₀::Vector{arrayType};
                 parameters::pType=nothing, periodicity=zero(q₀[begin])) where {
                        dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
                        vType <: Function, BType <: Function,
                        pType <: Union{NamedTuple,Nothing}}

        d  = length(q₀[begin])
        ni = length(q₀)

        @assert all(length(q) == d for q in q₀)

        @assert ns ≥ 1
        @assert ni == 1 || ns == 1
        # either multiple deterministic initial conditions and one sample path
        # or one deterministic initial condition and multiple sample paths

        new{dType,tType,arrayType,vType,BType,pType}(d, m, ns, v, B, t₀, q₀, parameters, periodicity)
    end
end


SDE(m, ns, v, B, q₀::StateVector; kwargs...) = SDE(m, ns, v, B, 0.0, q₀; kwargs...)
SDE(m, ns, v, B, t₀, q₀::State; kwargs...) = SDE(m, ns, v, B, t₀, [q₀]; kwargs...)
SDE(m, ns, v, B, q₀::State; kwargs...) = SDE(m, ns, v, B, 0.0, q₀; kwargs...)

# const SDEHT{HT,DT,TT,AT,VT,BT,PT} = SDE{DT,TT,AT,VT,BT,HT,PT} # type alias for dispatch on Hamiltonian type parameter
# const SDEPT{PT,DT,TT,AT,VT,BT,HT} = SDE{DT,TT,AT,VT,BT,HT,PT} # type alias for dispatch on parameters type parameter
const SDEPT{PT,DT,TT,AT,VT,BT} = SDE{DT,TT,AT,VT,BT,PT} # type alias for dispatch on parameters type parameter

Base.hash(sde::SDE, h::UInt) = hash(sde.d, hash(sde.m, hash(sde.ns,
                               hash(sde.v, hash(sde.B, hash(sde.t₀, hash(sde.q₀,
                               hash(sde.parameters, hash(sde.periodicity, h)))))))))

Base.:(==)(sde1::SDE, sde2::SDE) = (
                                sde1.d == sde2.d
                             && sde1.m == sde2.m
                             && sde1.ns == sde2.ns
                             && sde1.v == sde2.v
                             && sde1.B == sde2.B
                             && sde1.t₀ == sde2.t₀
                             && sde1.q₀ == sde2.q₀
                             && sde1.parameters == sde2.parameters
                             && sde1.periodicity == sde2.periodicity)

function Base.similar(equ::SDE, t₀::Real, q₀::StateVector, ns::Int=equ.ns)
    @assert all([length(q) == ndims(equ) for q in q₀])
    SDE(equ.m, ns, equ.v, equ.B, t₀, q₀; parameters=equ.parameters, periodicity=equ.periodicity)
end

Base.similar(equ::SDE, t₀::Real, q₀::State, ns::Int=equ.ns; kwargs...) = similar(equ, t₀, [q₀], ns; kwargs...)
Base.similar(equ::SDE, q₀::Union{State,StateVector}, ns::Int=equ.ns; kwargs...) = similar(equ, equ.t₀, q₀, ns; kwargs...)

Base.ndims(sde::SDE) = sde.d
Base.axes(equ::SDE) = axes(equ.q₀[begin])
Common.nsamples(equ::SDE) = length(equ.q₀)
Common.periodicity(equ::SDE) = equ.periodicity

initial_conditions(equ::SDE) = (equ.t₀, equ.q₀)

# hashamiltonian(::SDEHT{<:Nothing}) = false
# hashamiltonian(::SDEHT{<:Function}) = true

hasparameters(::SDEPT{<:Nothing}) = false
hasparameters(::SDEPT{<:NamedTuple}) = true

_get_v(equ::SDE) = hasparameters(equ) ? (t,q,v) -> equ.v(t, q, v, equ.parameters) : equ.v
_get_B(equ::SDE) = hasparameters(equ) ? (t,q,B,col=0) -> equ.B(t, q, B, equ.parameters, col) : equ.B
# _get_h(equ::SDE) = hasparameters(equ) ? (t,q) -> equ.h(t, q, equ.parameters) : equ.h

function get_function_tuple(equ::SDE)
    names = (:v,:B)
    equs  = (_get_v(equ), _get_B(equ))

    # if hashamiltonian(equ)
    #     names = (names..., :h)
    #     equs  = (equs..., _get_h(equ))
    # end

    NamedTuple{names}(equs)
end
