@doc raw"""
`DAE`: Differential Algebraic Equation

Defines a differential algebraic initial value problem
```math
\begin{align*}
\dot{q} (t) &= v(t, q(t)) + u(t, q(t), \lambda(t)) , & q(t_{0}) &= q_{0} , \\
0 &= \phi (t, q(t), \lambda(t)) , & \lambda(t_{0}) &= \lambda_{0} ,
\end{align*}
```
with vector field ``v``, projection ``u``, algebraic constraint ``\phi=0``,
initial conditions ``q_{0}`` and ``\lambda_{0}``, the dynamical variable ``q``
taking values in ``\mathbb{R}^{m}`` and the algebraic variable ``\lambda``
taking values in ``\mathbb{R}^{n}``.

### Fields

* `d`: dimension of dynamical variable ``q`` and the vector field ``v``
* `m`: dimension of algebraic variable ``\lambda`` and the constraint ``\phi``
* `n`: number of initial conditions
* `v`: function computing the vector field
* `u`: function computing the projection
* `ϕ`: algebraic constraint
* `t₀`: initial time
* `q₀`: initial condition for dynamical variable ``q``
* `λ₀`: initial condition for algebraic variable ``\lambda``


The function `v`, providing the vector field, takes three arguments,
`v(t, q, v)`, the functions `u` and `ϕ`, providing the projection and the
algebraic constraint take four arguments, `u(t, q, λ, u)` and `ϕ(t, q, λ, ϕ)`,
where `t` is the current time, `q` and `λ` are the current solution vectors,
and `v`, `u` and `ϕ` are the vectors which hold the result of evaluating the
vector field ``v``, the projection ``u`` and the algebraic constraint ``\phi``
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
struct DAE{dType <: Number, tType <: Number, vType <: Function, uType <: Function,
           ϕType <: Function, pType <: Union{Tuple,Nothing}, N} <: Equation{dType, tType}

    d::Int
    m::Int
    n::Int
    v::vType
    u::uType
    ϕ::ϕType
    t₀::tType
    q₀::Array{dType, N}
    λ₀::Array{dType, N}
    parameters::pType
    periodicity::Vector{dType}

    function DAE(DT::DataType, N::Int, d::Int, m::Int, n::Int,
                 v::vType, u::uType, ϕ::ϕType, t₀::tType,
                 q₀::DenseArray{dType}, λ₀::DenseArray{dType};
                 parameters=nothing, periodicity=zeros(DT,d)) where {
                        dType <: Number, tType <: Number, vType <: Function,
                        uType <: Function, ϕType <: Function}

        @assert d == size(q₀,1)
        @assert m == size(λ₀,1)
        @assert n == size(q₀,2) == size(λ₀,2)
        @assert d ≥ m
        @assert ndims(q₀) == ndims(λ₀) == N ∈ (1,2)

        new{DT, tType, vType, uType, ϕType, typeof(parameters), N}(d, m, n, v, u, ϕ, t₀,
                convert(Array{DT}, q₀), convert(Array{DT}, λ₀), parameters, periodicity)
    end
end

function DAE(v, u, ϕ, t₀, q₀::DenseArray{DT}, λ₀::DenseArray{DT}; kwargs...) where {DT}
    DAE(DT, ndims(q₀), size(q₀,1), size(λ₀,1), size(q₀,2), v, u, ϕ, t₀, q₀, λ₀; kwargs...)
end

function DAE(v, u, ϕ, q₀, λ₀; kwargs...)
    DAE(v, u, ϕ, zero(eltype(q₀)), q₀, λ₀; kwargs...)
end

Base.hash(dae::DAE, h::UInt) = hash(dae.d, hash(dae.m, hash(dae.n, hash(dae.v,
        hash(dae.u, hash(dae.t₀, hash(dae.q₀, hash(dae.λ₀,
        hash(dae.periodicity, hash(dae.parameters, h))))))))))

Base.:(==)(dae1::DAE, dae2::DAE) = (
                                dae1.d == dae2.d
                             && dae1.m == dae2.m
                             && dae1.n == dae2.n
                             && dae1.v == dae2.v
                             && dae1.u == dae2.u
                             && dae1.t₀ == dae2.t₀
                             && dae1.q₀ == dae2.q₀
                             && dae1.λ₀ == dae2.λ₀
                             && dae1.parameters == dae1.parameters
                             && dae1.periodicity == dae1.periodicity)

function Base.similar(dae::DAE, q₀, λ₀=get_λ₀(q₀, dae.λ₀); kwargs...)
    similar(dae, dae.t₀, q₀, λ₀; kwargs...)
end

function Base.similar(dae::DAE, t₀::TT, q₀::DenseArray{DT}, λ₀::DenseArray{DT}=get_λ₀(q₀, dae.λ₀);
                      parameters=dae.parameters, periodicity=dae.periodicity) where {DT  <: Number, TT <: Number}
    @assert dae.d == size(q₀,1)
    @assert dae.m == size(λ₀,1)
    DAE(dae.v, dae.u, dae.ϕ, t₀, q₀, λ₀; parameters=parameters, periodicity=periodicity)
end

Base.ndims(dae::DAE) = dae.d
