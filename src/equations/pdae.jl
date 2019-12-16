@doc raw"""
`PDAE`: Partitioned Differential Algebraic Equation

Defines a partitioned differential algebraic initial value problem
```math
\begin{align*}
\dot{q} (t) &= v(t, q(t), p(t)) + u(t, q(t), p(t), \lambda(t)) , & q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), p(t)) + r(t, q(t), p(t), \lambda(t)) , & p(t_{0}) &= p_{0} , \\
0 &= \phi (t, q(t), p(t), \lambda(t)) , & \lambda(t_{0}) &= \lambda_{0} ,
\end{align*}
```
with vector fields ``v`` and ``f``, projection ``u`` and ``r``,
algebraic constraint ``\phi=0``,
conditions ``(q_{0}, p_{0})`` and ``\lambda_{0}``, the dynamical variables
``(q,p)`` taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}`` and
the algebraic variable ``\lambda`` taking values in ``\mathbb{R}^{n}``.

### Fields

* `d`: dimension of dynamical variables ``q`` and ``p`` as well as the vector fields ``f`` and ``p``
* `m`: dimension of algebraic variable ``\lambda`` and the constraint ``\phi``
* `n`: number of initial conditions
* `v`: function computing the vector field ``v``
* `f`: function computing the vector field ``f``
* `u`: function computing the projection
* `g`: function computing the projection
* `ϕ`: algebraic constraint
* `t₀`: initial time
* `q₀`: initial condition for dynamical variable ``q``
* `p₀`: initial condition for dynamical variable ``p``
* `λ₀`: initial condition for algebraic variable ``\lambda``

"""
struct PDAE{dType <: Number, tType <: Number,
            vType <: Function, fType <: Function,
            uType <: Function, gType <: Function,
            ϕType <: Function, pType <: Union{Tuple,Nothing}, N} <: AbstractEquationPDAE{dType, tType}

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
    parameters::pType
    periodicity::Vector{dType}

    function PDAE(DT::DataType, N::Int, d::Int, m::Int, n::Int,
                  v::vType, f::fType, u::uType, g::gType, ϕ::ϕType,
                  t₀::tType, q₀::DenseArray{dType}, p₀::DenseArray{dType}, λ₀::DenseArray{dType};
                  parameters=nothing, periodicity=zeros(DT,d)) where {
                         dType <: Number, tType <: Number,
                         vType <: Function, fType <: Function,
                         uType <: Function, gType <: Function,
                         ϕType <: Function}

        @assert d == size(q₀,1) == size(p₀,1)
        @assert m == size(λ₀,1)
        @assert n == size(q₀,2) == size(p₀,2) == size(λ₀,2)
        @assert 2n ≥ m
        @assert ndims(q₀) == ndims(p₀) == ndims(λ₀) == N ∈ (1,2)

        new{DT, tType, vType, fType, uType, gType, ϕType, typeof(parameters), N}(d, m, n, v, f, u, g, ϕ, t₀,
                convert(Array{DT}, q₀), convert(Array{DT}, p₀), convert(Array{DT}, λ₀),
                parameters, periodicity)
    end
end

function PDAE(v, f, u, g, ϕ, t₀, q₀::DenseArray{DT}, p₀::DenseArray{DT}, λ₀::DenseArray{DT}; kwargs...) where {DT}
    PDAE(DT, ndims(q₀), size(q₀,1), size(λ₀,1), size(q₀,2), v, f, u, g, ϕ, t₀, q₀, p₀, λ₀; kwargs...)
end

function PDAE(v, f, u, g, ϕ, q₀, p₀, λ₀; kwargs...)
    PDAE(v, f, u, g, ϕ, zero(eltype(q₀)), q₀, p₀, λ₀; kwargs...)
end

Base.hash(dae::PDAE, h::UInt) = hash(dae.d, hash(dae.m, hash(dae.n,
        hash(dae.v, hash(dae.f, hash(dae.u, hash(dae.g, hash(dae.t₀,
        hash(dae.q₀, hash(dae.p₀, hash(dae.λ₀, hash(dae.periodicity, hash(dae.parameters, h)))))))))))))

Base.:(==)(dae1::PDAE, dae2::PDAE) = (
                                dae1.d == dae2.d
                             && dae1.m == dae2.m
                             && dae1.n == dae2.n
                             && dae1.v == dae2.v
                             && dae1.f == dae2.f
                             && dae1.u == dae2.u
                             && dae1.g == dae2.g
                             && dae1.ϕ == dae2.ϕ
                             && dae1.t₀ == dae2.t₀
                             && dae1.q₀ == dae2.q₀
                             && dae1.p₀ == dae2.p₀
                             && dae1.λ₀ == dae2.λ₀
                             && dae1.parameters == dae1.parameters
                             && dae1.periodicity == dae1.periodicity)

function Base.similar(dae::PDAE, q₀, p₀, λ₀=get_λ₀(q₀, dae.λ₀); kwargs...)
    similar(dae, dae.t₀, q₀, p₀, λ₀; kwargs...)
end

function Base.similar(dae::PDAE, t₀::TT, q₀::DenseArray{DT}, p₀::DenseArray{DT}, λ₀::DenseArray{DT}=get_λ₀(q₀, dae.λ₀);
                      parameters=dae.parameters, periodicity=dae.periodicity) where {DT  <: Number, TT <: Number}
    @assert dae.d == size(q₀,1) == size(p₀,1)
    @assert dae.m == size(λ₀,1)
    PDAE(dae.v, dae.f, dae.u, dae.g, dae.ϕ, t₀, q₀, p₀, λ₀; parameters=parameters, periodicity=periodicity)
end

@inline Base.ndims(dae::PDAE) = dae.d

@inline periodicity(equation::PDAE) = equation.periodicity
