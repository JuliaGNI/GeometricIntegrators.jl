
"""
`SDE`: Stochastic Differential Equation

Defines a stochastic differential initial value problem
```math
\\begin{align*}
\\dq (t) &= v(t, q(t)) \, dt + u(t, q(t)), dW , & q(t_{0}) &= q_{0} ,
\\end{align*}
```
with deterministic vector field ``v``, stochastic vector field ``u``,
initial conditions ``q_{0}`` and the dynamical variable ``q``
taking values in ``\\mathbb{R}^{m}``.

### Fields

* `d`: dimension of dynamical variable ``q`` and the vector field ``v``
* `n`: number of initial conditions
* `v`: function computing the deterministic vector field
* `u`: function computing the stochastic vector field
* `t₀`: initial time
* `q₀`: initial condition for dynamical variable ``q``


The functions `v` and `u`, providing the vector fields, take three arguments,
`v(t, q, v)` and `u(t, q, u)`, where `t` is the current time, `q` is the
current solution vector, and `v` and `u` are the vectors which hold the result
of evaluating the vector fields ``v`` and ``u`` on `t` and `q`.

### Example

```julia
    function v(λ, t, q, v)
        v[1] = λ*q[1]
        v[2] = λ*q[2]
    end

    function u(μ, t, q, u)
        u[1] = μ*q[1]
        u[2] = μ*q[2]
    end

    t₀ = 0.
    q₀ = [1., 1.]
    λ  = 2.
    μ  = 1.

    v_sde = (t, q, v) -> v(λ, t, q, v)
    u_sde = (t, q, v) -> u(μ, t, q, v)

    sde = SDE(v_sde, u_sde, t₀, q₀)
```
"""
struct SDE{dType <: Number, tType <: Number, vType <: Function, uType <: Function, N} <: Equation{dType, tType}
    d::Int
    n::Int
    v::vType
    u::uType
    t₀::tType
    q₀::Array{dType, N}

    function SDE{dType,tType,vType,uType,N}(d, n, v, u, t₀, q₀) where {dType <: Number, tType <: Number, vType <: Function, uType <: Function, N}
        @assert d == size(q₀,1)
        @assert n == size(q₀,2)

        @assert dType == eltype(q₀)
        @assert tType == typeof(t₀)

        @assert ndims(q₀) == N ∈ (1,2)

        new(d, n, v, u, t₀, q₀)
    end
end

function SDE(v::VT, u::UT, t₀::TT, q₀::DenseArray{DT}) where {DT,TT,VT,UT}
    SDE{DT, TT, VT, UT, ndims(q₀)}(size(q₀, 1), size(q₀, 2), v, u, t₀, q₀)
end

function SDE(v, u, q₀)
    SDE(v, u, zero(eltype(q₀)), q₀)
end

Base.hash(sde::SDE, h::UInt) = hash(sde.d, hash(sde.n, hash(sde.v, hash(sde.u, hash(sde.t₀, hash(sde.q₀, h))))))
Base.:(==)(sde1::SDE, sde2::SDE) = (
                                sde1.d == sde2.d
                             && sde1.n == sde2.n
                             && sde1.v == sde2.v
                             && sde1.u == sde2.u
                             && sde1.t₀ == sde2.t₀
                             && sde1.q₀ == sde2.q₀)

function Base.similar(sde::SDE{DT,TT,VT,UT}, q₀::DenseArray{DT}) where {DT, TT, VT, UT}
    similar(sde, sde.t₀, q₀)
end

function Base.similar(sde::SDE{DT,TT,VT,UT}, t₀::TT, q₀::DenseArray{DT}) where {DT, TT, VT, UT}
    @assert sde.d == size(q₀,1)
    SDE(sde.v, t₀, q₀)
end
