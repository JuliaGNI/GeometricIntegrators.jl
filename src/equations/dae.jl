
"""
`DAE`: Differential Algebraic Equation

Defines a differential algebraic initial value problem
```math
\\begin{align*}
\\dot{q} (t) &= f(t, q(t)) + u(t, q(t), \\lambda(t)) , & q(t_{0}) &= q_{0} , \\\\
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
* `f`: function computing the vector field
* `u`: function computing the projection
* `ϕ`: algebraic constraint
* `t₀`: initial time
* `q₀`: initial condition for dynamical variable ``q``
* `λ₀`: initial condition for algebraic variable ``\\lambda``


The function `f`, providing the vector field, takes three arguments,
`f(t, q, fq)`, the functions `u` and `ϕ`, providing the projection and the
algebraic constraint take four arguments, `u(t, q, λ, fu)` and `ϕ(t, q, λ, fϕ)`,
where `t` is the current time, `q` and `λ` are the current solution vectors,
and `fq`, `fu` and `fϕ` are the vectors which hold the result of evaluating the
vector field ``f``, the projection ``u`` and the algebraic constraint ``\\phi``
on `t`, `q` and `λ`.

### Example

```julia
    function f(t, q, fq)
        fq[1] = q[1]
        fq[2] = q[2]
    end

    function u(t, q, λ, fu)
        fu[1] = +λ[1]
        fu[2] = -λ[1]
    end

    function ϕ(t, q, λ, fϕ)
        fϕ[1] = q[2] - q[1]
    end

    t₀ = 0.
    q₀ = [1., 1.]
    λ₀ = [0.]

    dae = DAE(f, u, ϕ, t₀, q₀, λ₀)

```
"""
immutable DAE{T} <: Equation{T}
    m::Int
    n::Int
    f::Function
    u::Function
    ϕ::Function
    t₀::T
    q₀::Array{T,1}
    λ₀::Array{T,1}

    function DAE(m, n, f, u, ϕ, t₀, q₀, λ₀)
        @assert m == length(q₀)
        @assert n == length(λ₀)
        @assert m ≥ n
        @assert T == eltype(q₀)
        @assert T == eltype(λ₀)

        new(m, n, f, u, ϕ, t₀, q₀, λ₀)
    end
end

function DAE{T}(f::Function, u::Function, ϕ::Function, t₀::Real, q₀::Vector{T}, λ₀::Vector{T})
    DAE{T}(length(q₀), length(λ₀), f, u, ϕ, t₀, q₀, λ₀)
end

function DAE{T}(f::Function, u::Function, ϕ::Function, q₀::Vector{T}, λ₀::Vector{T})
    DAE(f, u, ϕ, 0, q₀, λ₀)
end
