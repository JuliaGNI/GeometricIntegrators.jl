
"""
`DAE`: Differential Algebraic Equation

Defines a differential algebraic initial value problem
```math
\\begin{align*}
\\dot{q} (t) &= v(t, q(t)) + u(t, q(t), \\lambda(t)) , & q(t_{0}) &= q_{0} , \\\\
0 &= \\phi (t, q(t), \\lambda(t)) , & \\lambda(t_{0}) &= \\lambda_{0} ,
\\end{align*}
```
with vector field ``f``, projection ``u``, algebraic constraint ``\\phi=0``,
initial conditions ``q_{0}`` and ``\\lambda_{0}``, the dynamical variable ``q``
taking values in ``\\mathbb{R}^{m}`` and the algebraic variable ``\\lambda``
taking values in ``\\mathbb{R}^{n}``.

### Fields

* `m`: dimension of dynamical variable ``q`` and the vector field ``f``
* `n`: dimension of algebraic variable ``\\lambda`` and the constraint ``\\phi``
* `v`: function computing the vector field
* `u`: function computing the projection
* `ϕ`: algebraic constraint
* `t₀`: initial time
* `q₀`: initial condition for dynamical variable ``q``
* `λ₀`: initial condition for algebraic variable ``\\lambda``


The function `v`, providing the vector field, takes three arguments,
`v(t, q, v)`, the functions `u` and `ϕ`, providing the projection and the
algebraic constraint take four arguments, `u(t, q, λ, u)` and `ϕ(t, q, λ, ϕ)`,
where `t` is the current time, `q` and `λ` are the current solution vectors,
and `v`, `u` and `ϕ` are the vectors which hold the result of evaluating the
vector field ``v``, the projection ``u`` and the algebraic constraint ``\\phi``
on `t`, `q` and `λ`.

### Example

```julia
    function v(t, q, v)
        v[1] = q[1]
        v[2] = q[2]
    end

    function u(t, q, λ, u)
        u[1] = +λ[1]
        u[2] = -λ[1]
    end

    function ϕ(t, q, λ, ϕ)
        ϕ[1] = q[2] - q[1]
    end

    t₀ = 0.
    q₀ = [1., 1.]
    λ₀ = [0.]

    dae = DAE(v, u, ϕ, t₀, q₀, λ₀)

```
"""
immutable DAE{dType <: Number, tType <: Number, vType <: Function, uType <: Function, ϕType <: Function} <: Equation{dType, tType}
    d::Int
    m::Int
    n::Int
    v::vType
    u::uType
    ϕ::ϕType
    t₀::tType
    q₀::Array{dType, 2}
    λ₀::Array{dType, 2}

    function DAE(d, m, n, v, u, ϕ, t₀, q₀, λ₀)
        @assert d == size(q₀,1)
        @assert m == size(λ₀,1)
        @assert n == size(q₀,2) == size(λ₀,2)
        @assert d ≥ m

        @assert dType == eltype(q₀)
        @assert dType == eltype(λ₀)
        @assert tType == typeof(t₀)

        @assert ndims(q₀) ∈ (1,2)
        @assert ndims(λ₀) ∈ (1,2)

        if ndims(q₀) == 1
            q₀ = reshape(q₀, d, n)
        end

        if ndims(λ₀) == 1
            λ₀ = reshape(λ₀, m, n)
        end

        new(d, m, n, v, u, ϕ, t₀, q₀, λ₀)
    end
end

function DAE{DT, TT, VT, UT, ΦT}(v::VT, u::UT, ϕ::ΦT, t₀::TT, q₀::DenseArray{DT}, λ₀::DenseArray{DT})
    DAE{DT, TT, VT, UT, ΦT}(size(q₀, 1), size(q₀, 2), size(λ₀, 1), v, u, ϕ, t₀, q₀, λ₀)
end

function DAE(v, u, ϕ, q₀, λ₀)
    DAE(v, u, ϕ, zero(eltype(q₀)), q₀, λ₀)
end

Base.hash(dae::DAE, h::UInt) = hash(dae.d, hash(dae.m, hash(dae.n, hash(dae.v, hash(dae.u, hash(dae.t₀, hash(dae.q₀, hash(dae.λ₀, h))))))))
Base.:(==)(dae1::DAE, dae2::DAE) = (
                                dae1.d == dae2.d
                             && dae1.m == dae2.m
                             && dae1.n == dae2.n
                             && dae1.v == dae2.v
                             && dae1.u == dae2.u
                             && dae1.t₀ == dae2.t₀
                             && dae1.q₀ == dae2.q₀
                             && dae1.λ₀ == dae2.λ₀)
