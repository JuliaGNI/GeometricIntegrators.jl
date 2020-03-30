@doc raw"""
`SDE`: Stratonovich Stochastic Differential Equation

Defines a stochastic differential initial value problem
```math
\begin{align*}
\dq (t) &= v(t, q(t)) \, dt + B(t, q(t)) \circ dW , & q(t_{0}) &= q_{0} ,
\end{align*}
```
with drift vector field ``v``, diffusion matrix ``B``,
initial conditions ``q_{0}``, the dynamical variable ``q``
taking values in ``\mathbb{R}^{d}``, and the m-dimensional Wiener process W

### Fields

* `d`:  dimension of dynamical variable ``q`` and the vector field ``v``
* `m`:  dimension of the Wiener process
* `ni`: number of initial conditions
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
    function v(λ, t, q, v)
        v[1] = λ*q[1]
        v[2] = λ*q[2]
    end

    function B(μ, t, q, B; col=0)
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

    v_sde = (t, q, v) -> v(λ, t, q, v)
    B_sde = (t, q, B) -> B(μ, t, q, B)

    sde = SDE(v_sde, B_sde, t₀, q₀)
```
"""
struct SDE{dType <: Number, tType <: Number, vType <: Function, BType <: Function,
           pType <: Union{NamedTuple,Nothing}, N} <: AbstractEquationSDE{dType, tType}

    d::Int
    m::Int
    ni::Int
    ns::Int
    v::vType
    B::BType
    t₀::tType
    q₀::Array{dType, N}
    parameters::pType
    periodicity::Vector{dType}

    function SDE(m, ns, v::vType, B::BType, t₀::tType, q₀::AbstractArray{dType,N};
                 parameters::pType=nothing, periodicity=zeros(dType,size(q₀,1))) where {
                        dType <: Number, tType <: Number, vType <: Function, BType <: Function,
                        pType <: Union{NamedTuple,Nothing}, N}

        d  = size(q₀,1)
        ni = size(q₀,2)

        @assert N ∈ (1,2)
        @assert ni ≥ 1
        @assert ns ≥ 1
        @assert ni == 1 || ns == 1
        # either multiple deterministic initial conditions and one sample path
        # or one deterministic initial condition and multiple sample paths

        new{dType,tType,vType,BType,pType,N}(d, m, ni, ns, v, B, t₀, q₀, parameters, periodicity)
    end
end


function SDE(m::Int, ns::Int, v::Function, B::Function, q₀::AbstractArray{DT}; kwargs...) where {DT}
    SDE(m, ns, v, B, zero(DT), q₀; kwargs...)
end


Base.hash(sde::SDE, h::UInt) = hash(sde.d, hash(sde.m, hash(sde.ni, hash(sde.ns,
                               hash(sde.v, hash(sde.B, hash(sde.t₀, hash(sde.q₀,
                               hash(sde.parameters, hash(sde.periodicity, h))))))))))

Base.:(==)(sde1::SDE, sde2::SDE) = (
                                sde1.d == sde2.d
                             && sde1.m == sde2.m
                             && sde1.ni == sde2.ni
                             && sde1.ns == sde2.ns
                             && sde1.v == sde2.v
                             && sde1.B == sde2.B
                             && sde1.t₀ == sde2.t₀
                             && sde1.q₀ == sde2.q₀
                             && sde1.parameters == sde2.parameters
                             && sde1.periodicity == sde2.periodicity)

function Base.similar(sde::SDE, q₀::AbstractArray)
    similar(sde, sde.t₀, q₀)
end

function Base.similar(sde::SDE, q₀::AbstractArray, ns::Int)
    similar(sde, sde.t₀, q₀, ns)
end

function Base.similar(sde::SDE, t₀::TT, q₀::AbstractArray{DT,N}, ns::Int=(N > 1 ? 1 : sde.ns)) where {DT <: Number, TT <: Number, N}
    @assert sde.d == size(q₀,1)
    SDE(sde.m, ns, sde.v, sde.B, t₀, q₀; parameters=sde.parameters, periodicity=sde.periodicity)
end

@inline Base.ndims(sde::SDE) = sde.d

@inline CommonFunctions.periodicity(equation::SDE) = equation.periodicity

function get_function_tuple(equation::SDE{DT,TT,VT,BT,Nothing}) where {DT, TT, VT, BT}
    NamedTuple{(:v,:B)}((equation.v, equation.B))
end

function get_function_tuple(equation::SDE{DT,TT,VT,BT,PT}) where {DT, TT, VT, BT, PT <: NamedTuple}
    vₚ = (t,q,v) -> equation.v(t, q, v, equation.parameters)
    Bₚ = (t,q,B,col=0) -> equation.B(t, q, B, equation.parameters, col)

    names = (:v,:B)
    equs  = (vₚ,Bₚ)

    NamedTuple{names}(equs)
end
