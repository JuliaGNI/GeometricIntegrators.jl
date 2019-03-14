@doc raw"""
`IDAE`: Implicit Differential Algebraic Equation

Defines a partitioned differential algebraic initial value problem
```math
\begin{align*}
\dot{q} (t) &= v(t) + u(t, q(t), p(t), \lambda(t)) , & q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), v(t)) + r(t, q(t), p(t), \lambda(t)) , & p(t_{0}) &= p_{0} , \\
p(t) &= p(t, q(t), v(t)) , && \
0 &= \phi (t, q(t), p(t), \lambda(t)) , & \lambda(t_{0}) &= \lambda_{0} ,
\end{align*}
```
with vector field ``f``, the momentum defined by ``p``, projection ``u`` and ``r``,
algebraic constraint ``\phi=0``,
conditions ``(q_{0}, p_{0})`` and ``\lambda_{0}``, the dynamical variables
``(q,p)`` taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}`` and
the algebraic variable ``\lambda`` taking values in ``\mathbb{R}^{n}``.

### Fields

* `d`: dimension of dynamical variables ``q`` and ``p`` as well as the vector fields ``f`` and ``p``
* `m`: dimension of algebraic variable ``\lambda`` and the constraint ``\phi``
* `n`: number of initial conditions
* `f`: function computing the vector field ``f``
* `p`: function computing ``p``
* `u`: function computing the projection
* `g`: function computing the projection
* `ϕ`: algebraic constraint
* `t₀`: initial time
* `q₀`: initial condition for dynamical variable ``q``
* `p₀`: initial condition for dynamical variable ``p``
* `λ₀`: initial condition for algebraic variable ``\lambda``

"""
struct IDAE{dType <: Number, tType <: Number, fType <: Function, pType <: Function, uType <: Function, gType <: Function, ϕType <: Function, vType <: Function, N} <: Equation{dType, tType}
    d::Int
    m::Int
    n::Int
    f::fType
    p::pType
    u::uType
    g::gType
    ϕ::ϕType
    v::vType
    t₀::tType
    q₀::Array{dType, N}
    p₀::Array{dType, N}
    λ₀::Array{dType, N}

    function IDAE{dType,tType,fType,pType,uType,gType,ϕType,vType,N}(d, m, n, f, p, u, g, ϕ, v, t₀, q₀, p₀, λ₀) where {dType <: Number, tType <: Number, fType <: Function, pType <: Function, uType <: Function, gType <: Function, ϕType <: Function, vType <: Function, N}
        @assert d == size(q₀,1) == size(p₀,1)
        @assert m == size(λ₀,1)
        @assert n == size(q₀,2) == size(p₀,2) == size(λ₀,2)
        @assert 2n ≥ m

        @assert dType == eltype(q₀)
        @assert dType == eltype(p₀)
        @assert dType == eltype(λ₀)

        @assert ndims(q₀) == ndims(p₀) == ndims(λ₀) == N ∈ (1,2)

        new(d, m, n, f, p, u, g, ϕ, v, t₀, q₀, p₀, λ₀)
    end
end

function IDAE(f::FT, p::PT, u::UT, g::GT, ϕ::ΦT, v::VT, t₀::TT, q₀::DenseArray{DT}, p₀::DenseArray{DT}, λ₀::DenseArray{DT}) where {DT,TT,FT,PT,UT,GT,ΦT,VT}
    @assert size(q₀) == size(p₀)
    @assert size(q₀,2) == size(λ₀,2)
    IDAE{DT, TT, FT, PT, UT, GT, ΦT, VT, ndims(q₀)}(size(q₀, 1), size(λ₀, 1), size(q₀, 2), f, p, u, g, ϕ, v, t₀, q₀, p₀, λ₀)
end

function IDAE(f::Function, p::Function, u::Function, g::Function, ϕ::Function, t₀::Number, q₀, p₀, λ₀)
    IDAE(f, p, u, g, ϕ, function_v_dummy, t₀, q₀, p₀, λ₀)
end

function IDAE(f::Function, p::Function, u::Function, g::Function, ϕ::Function, v::Function, q₀, p₀, λ₀)
    IDAE(f, p, u, g, ϕ, v, zero(Float64), q₀, p₀, λ₀)
end

function IDAE(f, p, u, g, ϕ, q₀, p₀, λ₀)
    IDAE(f, p, u, g, ϕ, function_v_dummy, zero(Float64), q₀, p₀, λ₀)
end

Base.hash(dae::IDAE, h::UInt) = hash(dae.d, hash(dae.m, hash(dae.n, hash(dae.f, hash(dae.p, hash(dae.u, hash(dae.g, hash(dae.ϕ, hash(dae.v, hash(dae.t₀, hash(dae.q₀, hash(dae.p₀, hash(dae.λ₀, h)))))))))))))
Base.:(==)(dae1::IDAE, dae2::IDAE) = (
                                dae1.d == dae2.d
                             && dae1.m == dae2.m
                             && dae1.n == dae2.n
                             && dae1.f == dae2.f
                             && dae1.p == dae2.p
                             && dae1.u == dae2.u
                             && dae1.g == dae2.g
                             && dae1.ϕ == dae2.ϕ
                             && dae1.v == dae2.v
                             && dae1.t₀ == dae2.t₀
                             && dae1.q₀ == dae2.q₀
                             && dae1.p₀ == dae2.p₀
                             && dae1.λ₀ == dae2.λ₀)
