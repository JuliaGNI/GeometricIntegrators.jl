@doc raw"""
`SPDAE`: Split Partitioned Differential Algebraic Equation *EXPERIMENTAL*

Defines a split differential algebraic initial value problem, that is
a canonical Hamiltonian system of equations subject to Dirac constraints,
```math
\begin{aligned}
\dot{q} (t) &= v_1(t, q(t), p(t)) + v_2(t, q(t), p(t), \lambda(t)) + v_3(t, q(t), p(t), \lambda(t), \gamma(t)) , & q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f_1(t, q(t), p(t)) + f_2(t, q(t), p(t), \lambda(t)) + f_3(t, q(t), p(t), \lambda(t), \gamma(t)) , & p(t_{0}) &= p_{0} , \\
0 &= \phi (t, q(t), p(t)) , \\
0 &= \psi (t, q(t), p(t), \dot{q}(t), \dot{p}(t)) ,
\end{aligned}
```
with vector fields ``v_i`` and ``f_i`` for ``i = 1 ... 3``,
primary constraint ``\phi(q,p)=0`` and secondary constraint ``\psi(q,p,\lambda)=0``,
initial conditions ``(q_{0}, p_{0})``, the dynamical variables ``(q,p)``
taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}`` and
the algebraic variables ``(\lambda, \gamma)`` taking values in
``\mathbb{R}^{n} \times \mathbb{R}^{d}``.

### Fields

* `d`: dimension of dynamical variables ``q`` and ``p`` as well as the vector fields ``v`` and ``f``
* `m`: dimension of algebraic variables ``\lambda`` and ``\gamma`` and the constraints ``\phi`` and ``\psi``
* `n`: number of initial conditions
* `v`: tuple of functions computing the vector fields ``v_i``, ``i = 1 ... 3``
* `f`: tuple of functions computing the vector fields ``f_i``, ``i = 1 ... 3``
* `ϕ`: primary constraints
* `ψ`: secondary constraints
* `t₀`: initial time
* `q₀`: initial condition for dynamical variable ``q``
* `p₀`: initial condition for dynamical variable ``p``
* `λ₀`: initial condition for algebraic variable ``λ``
"""
struct SPDAE{dType <: Number, tType <: Number, vType <: Tuple, fType <: Tuple,
             ϕType <: Function, ψType <: Function, pType <: Union{NamedTuple,Nothing},
             N} <: AbstractEquationPDAE{dType, tType}
    d::Int
    m::Int
    n::Int
    v::vType
    f::fType
    ϕ::ϕType
    ψ::ψType
    t₀::tType
    q₀::Array{dType, N}
    p₀::Array{dType, N}
    λ₀::Array{dType, N}
    parameters::pType
    periodicity::Vector{dType}

    function SPDAE(DT::DataType, N::Int, d::Int, m::Int, n::Int,
                   v::vType, f::fType, ϕ::ϕType, ψ::ψType, t₀::tType,
                   q₀::AbstractArray{dType}, p₀::AbstractArray{dType}, λ₀::AbstractArray{dType};
                   parameters::pType=nothing, periodicity=zeros(DT,d)) where {
                        dType <: Number, tType <: Number,
                        vType <: Tuple, fType <: Tuple,
                        ϕType <: Function, ψType <: Function,
                        pType <: Union{NamedTuple,Nothing}}

        @assert d == size(q₀,1) == size(p₀,1)
        @assert n == size(q₀,2) == size(p₀,2)
        @assert ndims(q₀) == ndims(p₀) == N ∈ (1,2)
        @assert 2d ≥ m

        new{DT, tType, vType, fType, ϕType, ψType, pType, N}(d, m, n, v, f, ϕ, ψ, t₀,
                convert(Array{DT}, q₀), convert(Array{DT}, p₀), convert(Array{DT}, λ₀),
                parameters, periodicity)
    end
end

function SPDAE(v, f, ϕ, ψ, t₀::Number, q₀::AbstractArray, p₀::AbstractArray, λ₀::AbstractArray=zero(q₀); kwargs...)
    SPDAE(eltype(q₀), ndims(q₀), size(q₀,1), size(λ₀,1), size(q₀,2), v, f, ϕ, ψ, t₀, q₀, p₀, λ₀; kwargs...)
end

function SPDAE(v, f, ϕ, ψ, q₀::AbstractArray, p₀::AbstractArray, λ₀::AbstractArray=zero(q₀); kwargs...)
    SPDAE(v, f, ϕ, ψ, zero(eltype(q₀)), q₀, p₀, λ₀; kwargs...)
end

Base.hash(dae::SPDAE, h::UInt) = hash(dae.d, hash(dae.m, hash(dae.n,
                    hash(dae.v, hash(dae.f, hash(dae.ϕ, hash(dae.ψ,
                    hash(dae.t₀, hash(dae.q₀, hash(dae.p₀, h))))))))))

Base.:(==)(dae1::SPDAE, dae2::SPDAE) = (
                                dae1.d == dae2.d
                             && dae1.m == dae2.m
                             && dae1.n == dae2.n
                             && dae1.v == dae2.v
                             && dae1.f == dae2.f
                             && dae1.ϕ == dae2.ϕ
                             && dae1.ψ == dae2.ψ
                             && dae1.t₀ == dae2.t₀
                             && dae1.q₀ == dae2.q₀
                             && dae1.p₀ == dae2.p₀)

function Base.similar(dae::SPDAE, q₀, p₀, λ₀=get_λ₀(q₀, dae.λ₀); kwargs...)
    similar(dae, dae.t₀, q₀, p₀, λ₀; kwargs...)
end

function Base.similar(dae::SPDAE, t₀::TT, q₀::AbstractArray{DT}, p₀::AbstractArray{DT}, λ₀::AbstractArray{DT}=get_λ₀(q₀, dae.λ₀);
                      parameters=dae.parameters, periodicity=dae.periodicity) where {DT  <: Number, TT <: Number}
    @assert dae.d == size(q₀,1) == size(p₀,1)
    @assert dae.m == size(λ₀,1)
    SPDAE(dae.v, dae.f, dae.ϕ, dae.ψ, t₀, q₀, p₀, λ₀;
         parameters=parameters, periodicity=periodicity)
end

@inline Base.ndims(equation::SPDAE) = equation.d
@inline CommonFunctions.nconstraints(equation::SPDAE) = equation.m
@inline CommonFunctions.periodicity(equation::SPDAE) = equation.periodicity

function get_function_tuple(equation::SPDAE{DT,TT,VT,FT,ϕT,ψT,Nothing}) where {DT, TT, VT, FT, ϕT, ψT}
    NamedTuple{(:v, :f, :ϕ, :ψ)}((equation.v, equation.f, equation.ϕ, equation.ψ))
end

function get_function_tuple(equation::SPDAE{DT,TT,VT,FT,ϕT,ψT,PT}) where {DT, TT, VT, FT, ϕT, ψT, PT <: NamedTuple}
    vₚ = (t,q,p,v)   -> equation.v(t, q, p, v, equation.parameters)
    fₚ = (t,q,p,f)   -> equation.f(t, q, p, f, equation.parameters)
    ϕₚ = (t,q,p,ϕ)   -> equation.ϕ(t, q, p, ϕ, equation.parameters)
    ψₚ = (t,q,p,v,f,ψ) -> equation.ψ(t, q, p, v, f, ψ, equation.parameters)

    names = (:v, :f, :ϕ, :ψ)
    equs  = (vₚ, fₚ, ϕₚ, ψₚ)

    NamedTuple{names}(equs)
end
