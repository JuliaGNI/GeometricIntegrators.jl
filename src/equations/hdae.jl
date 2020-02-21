@doc raw"""
`HDAE`: Hamiltonian Differential Algebraic Equation *EXPERIMENTAL*

Defines a Hamiltonian differential algebraic initial value problem, that is
a canonical Hamiltonian system of equations subject to Dirac constraints,
```math
\begin{align*}
\dot{q} (t) &= v_1(t, q(t), p(t)) + v_2(t, q(t), p(t), \lambda(t)) + v_3(t, q(t), p(t), \lambda(t), \gamma(t)) , & q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f_1(t, q(t), p(t)) + f_2(t, q(t), p(t), \lambda(t)) + f_3(t, q(t), p(t), \lambda(t), \gamma(t)) , & p(t_{0}) &= p_{0} , \\
0 &= \phi (t, q(t), p(t)) , \\
0 &= \psi (t, q(t), p(t), \dot{q}(t), \dot{p}(t)) ,
\end{align*}
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

"""
struct HDAE{dType <: Number, tType <: Number, vType <: Tuple, fType <: Tuple, ϕType <: Function, ψType <: Function, N} <: AbstractEquationPDAE{dType, tType}
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

    function HDAE{dType,tType,vType,fType,ϕType,ψType,N}(d, m, n, v, f, ϕ, ψ, t₀, q₀, p₀) where {dType <: Number, tType <: Number, vType <: Tuple, fType <: Tuple, ϕType <: Function, ψType <: Function, N}
        @assert d == size(q₀,1) == size(p₀,1)
        @assert n == size(q₀,2) == size(p₀,2)
        @assert 2n ≥ m

        @assert dType == eltype(q₀)
        @assert dType == eltype(p₀)

        @assert ndims(q₀) == ndims(p₀) == N ∈ (1,2)

        new(d, m, n, v, f, ϕ, ψ, t₀, q₀, p₀)
    end
end

function HDAE(v::VT, f::FT, ϕ::ΦT, ψ::ΨT, m::Int, t₀::TT, q₀::AbstractArray{DT}, p₀::AbstractArray{DT}) where {DT,TT,VT,FT,ΦT,ΨT}
    @assert size(q₀) == size(p₀)
    HDAE{DT, TT, VT, FT, ΦT, ΨT, ndims(q₀)}(size(q₀, 1), m, size(q₀, 2), v, f, ϕ, ψ, t₀, q₀, p₀)
end

function HDAE(v::Function, f::Function, ϕ::Function, ψ::Function, t₀, q₀, p₀)
    HDAE(v, f, ϕ, ψ, size(q₀,1), t₀, q₀, p₀)
end

function HDAE(v::Function, f::Function, ϕ::Function, ψ::Function, q₀, p₀)
    HDAE(v, f, ϕ, ψ, size(q₀,1), zero(Float64), q₀, p₀)
end

Base.hash(dae::HDAE, h::UInt) = hash(dae.d, hash(dae.m, hash(dae.n, hash(dae.v, hash(dae.f, hash(dae.ϕ, hash(dae.ψ, hash(dae.t₀, hash(dae.q₀, hash(dae.p₀, h))))))))))
Base.:(==)(dae1::HDAE, dae2::HDAE) = (
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

@inline Base.ndims(dae::HDAE) = ode.d

@inline CommonFunctions.periodicity(equation::HDAE) = equation.periodicity
