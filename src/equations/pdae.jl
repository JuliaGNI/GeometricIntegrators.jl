
"""
`PDAE`: Partitioned Differential Algebraic Equation

Defines a partitioned differential algebraic initial value problem
```math
\\begin{align*}
\\dot{q} (t) &= v(t, q(t), p(t)) + u(t, q(t), p(t), \\lambda(t)) , & q(t_{0}) &= q_{0} , \\\\
\\dot{p} (t) &= f(t, q(t), p(t)) + r(t, q(t), p(t), \\lambda(t)) , & p(t_{0}) &= p_{0} , \\\\
0 &= \\phi (t, q(t), p(t), \\lambda(t)) , & \\lambda(t_{0}) &= \\lambda_{0} ,
\\end{align*}
```
with vector fields ``v`` and ``f``, projection ``u`` and ``r``,
algebraic constraint ``\\phi=0``,
conditions ``(q_{0}, p_{0})`` and ``\\lambda_{0}``, the dynamical variables
``(q,p)`` taking values in ``\\mathbb{R}^{d} \\times \\mathbb{R}^{d}`` and
the algebraic variable ``\\lambda`` taking values in ``\\mathbb{R}^{n}``.

### Fields

* `m`: dimension of dynamical variables ``q`` and ``p`` as well as the vector fields ``v`` and ``f``
* `n`: dimension of algebraic variable ``\\lambda`` and the constraint ``\\phi``
* `v`: function computing the vector field ``v``
* `f`: function computing the vector field ``f``
* `u`: function computing the projection
* `g`: function computing the projection
* `ϕ`: algebraic constraint
* `t₀`: initial time
* `q₀`: initial condition for dynamical variable ``q``
* `p₀`: initial condition for dynamical variable `p`
* `λ₀`: initial condition for algebraic variable ``\\lambda``

"""
immutable PDAE{dType <: Number, tType <: Number, vType <: Function, fType <: Function, uType <: Function, gType <: Function, ϕType <: Function, N} <: Equation{dType, tType}
    d::Int
    m::Int
    n::Int
    v::vType
    f::fType
    u::uType
    g::gType
    ϕ::ϕType
    t₀::tType
    q₀::Array{dType, N}
    p₀::Array{dType, N}
    λ₀::Array{dType, N}

    function PDAE(d, m, n, v, f, u, g, ϕ, t₀, q₀, p₀, λ₀)
        @assert d == size(q₀,1) == size(p₀,1)
        @assert m == size(λ₀,1)
        @assert n == size(q₀,2) == size(p₀,2) == size(λ₀,2)
        @assert 2n ≥ m

        @assert dType == eltype(q₀)
        @assert dType == eltype(p₀)
        @assert dType == eltype(λ₀)

        @assert ndims(q₀) == ndims(p₀) == ndims(λ₀) == N ∈ (1,2)

        new(d, m, n, v, f, u, g, ϕ, t₀, q₀, p₀, λ₀)
    end
end

function PDAE{DT, TT, VT, FT, UT, GT, ΦT}(v::VT, f::FT, u::UT, g::GT, ϕ::ΦT, t₀::TT, q₀::DenseArray{DT}, p₀::DenseArray{DT}, λ₀::DenseArray{DT})
    @assert size(q₀) == size(p₀)
    @assert size(q₀,2) == size(λ₀,2)
    PDAE{DT, TT, VT, FT, UT, GT, ΦT, ndims(q₀)}(size(q₀, 1), size(λ₀, 1), size(q₀, 2), v, f, u, g, ϕ, t₀, q₀, p₀, λ₀)
end

function PDAE(v, f, u, g, ϕ, q₀, p₀, λ₀)
    PDAE(v, f, u, g, ϕ, zero(eltype(q₀)), q₀, p₀, λ₀)
end

Base.hash(dae::PDAE, h::UInt) = hash(dae.d, hash(dae.m, hash(dae.n, hash(dae.v, hash(dae.u, hash(dae.f, hash(dae.g, hash(dae.t₀, hash(dae.q₀, hash(dae.p₀, hash(dae.λ₀, h)))))))))))
Base.:(==)(dae1::PDAE, dae2::PDAE) = (
                                dae1.d == dae2.d
                             && dae1.m == dae2.m
                             && dae1.n == dae2.n
                             && dae1.v == dae2.v
                             && dae1.u == dae2.u
                             && dae1.f == dae2.f
                             && dae1.g == dae2.g
                             && dae1.t₀ == dae2.t₀
                             && dae1.q₀ == dae2.q₀
                             && dae1.p₀ == dae2.p₀
                             && dae1.λ₀ == dae2.λ₀)
