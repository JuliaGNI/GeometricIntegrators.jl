@doc raw"""
`HDAE`: Hamiltonian Differential Algebraic Equation *EXPERIMENTAL*

Defines a Hamiltonian differential algebraic initial value problem, that is
a canonical Hamiltonian system of equations subject to Dirac constraints,
```math
\begin{align*}
\dot{q} (t) &= v(t, q(t), p(t)) + u(t, q(t), p(t), \lambda(t)) + \bar{g}(t, q(t), p(t), \lambda(t), \gamma(t)) , & q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), p(t)) + g(t, q(t), p(t), \lambda(t)) + \bar{f}(t, q(t), p(t), \lambda(t), \gamma(t)) , & p(t_{0}) &= p_{0} , \\
0 &= \phi (t, q(t), p(t)) , \\
0 &= \psi (t, q(t), p(t), \dot{q}(t), \dot{p}(t)) ,
\end{align*}
```
with vector fields ``v``, ``u``, ``\bar{u}`` and ``f``, ``g``, ``\bar{g}``,
primary constraint ``\phi(q,p)=0`` and secondary constraint ``\psi(q,p,\lambda)=0``,
initial conditions ``(q_{0}, p_{0})``, the dynamical variables ``(q,p)``
taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}`` and
the algebraic variables ``(\lambda, \gamma)`` taking values in
``\mathbb{R}^{n} \times \mathbb{R}^{d}``.

### Fields

* `d`: dimension of dynamical variables ``q`` and ``p`` as well as the vector fields ``v`` and ``f``
* `m`: dimension of algebraic variables ``\lambda`` and ``\gamma`` and the constraints ``\phi`` and ``\psi``
* `n`: number of initial conditions
* `v`: function computing the Hamiltonian vector field ``v``
* `f`: function computing the Hamiltonian vector field ``f``
* `u`: function computing the primary projection field ``u``
* `g`: function computing the primary projection field ``g``
* `u̅`: function computing the secondary projection field ``u̅``
* `g̅`: function computing the secondary projection field ``g̅``
* `ϕ`: primary constraints
* `ψ`: secondary constraints
* `v̄`: function computing an initial guess for the velocity field ``v``` (optional)
* `f̄`: function computing an initial guess for the force field ``f`` (optional)
* `h`: function computing the Hamiltonian ``H``
* `t₀`: initial time
* `q₀`: initial condition for dynamical variable ``q``
* `p₀`: initial condition for dynamical variable ``p``
* `λ₀`: initial condition for algebraic variable ``λ``
"""
struct HDAE{dType <: Number, tType <: Number, vType <: Function, fType <: Function,
            uType <: Function, gType <: Function, u̅Type <: Function, g̅Type <: Function,
            ϕType <: Function, ψType <: Function, hType <: Function,
            v̄Type <: Function, f̄Type <: Function, pType <: Union{NamedTuple,Nothing},
            N} <: AbstractEquationPDAE{dType, tType}
    d::Int
    m::Int
    n::Int
    v::vType
    f::fType
    u::uType
    g::gType
    u̅::u̅Type
    g̅::g̅Type
    ϕ::ϕType
    ψ::ψType
    h::hType
    v̄::v̄Type
    f̄::f̄Type
    t₀::tType
    q₀::Array{dType, N}
    p₀::Array{dType, N}
    λ₀::Array{dType, N}
    parameters::pType
    periodicity::Vector{dType}

    function HDAE(DT::DataType, N::Int, d::Int, m::Int, n::Int,
               v::vType, f::fType, u::uType, g::gType, u̅::u̅Type, g̅::g̅Type,
               ϕ::ϕType, ψ::ψType, h::hType, t₀::tType,
               q₀::AbstractArray{dType}, p₀::AbstractArray{dType}, λ₀::AbstractArray{dType};
               v̄::v̄Type=v, f̄::f̄Type=f, parameters::pType=nothing, periodicity=zeros(DT,d)) where {
                    dType <: Number, tType <: Number,
                    vType <: Function, fType <: Function,
                    uType <: Function, gType <: Function,
                    u̅Type <: Function, g̅Type <: Function,
                    ϕType <: Function, ψType <: Function,
                    hType <: Function,
                    v̄Type <: Function, f̄Type <: Function,
                    pType <: Union{NamedTuple,Nothing}}

        @assert d == size(q₀,1) == size(p₀,1)
        @assert n == size(q₀,2) == size(p₀,2)
        @assert ndims(q₀) == ndims(p₀) == N ∈ (1,2)
        @assert 2d ≥ m

        new{DT, tType, vType, fType, uType, gType, u̅Type, g̅Type, ϕType, ψType, hType, v̄Type, f̄Type, pType, N}(
                d, m, n, v, f, u, g, u̅, g̅, ϕ, ψ, h, v̄, f̄, t₀,
                convert(Array{DT}, q₀), convert(Array{DT}, p₀), convert(Array{DT}, λ₀),
                parameters, periodicity)
    end
end

function HDAE(v, f, u, g, u̅, g̅, ϕ, ψ, h, t₀::Number, q₀::AbstractArray, p₀::AbstractArray, λ₀::AbstractArray=zero(q₀); kwargs...)
    HDAE(eltype(q₀), ndims(q₀), size(q₀,1), size(λ₀,1), size(q₀,2), v, f, u, g, u̅, g̅, ϕ, ψ, h, t₀, q₀, p₀, λ₀; kwargs...)
end

function HDAE(v, f, u, g, u̅, g̅, ϕ, ψ, h, q₀::AbstractArray, p₀::AbstractArray, λ₀::AbstractArray=zero(q₀); kwargs...)
    HDAE(v, f, u, g, u̅, g̅, ϕ, ψ, h, zero(eltype(q₀)), q₀, p₀, λ₀; kwargs...)
end

Base.hash(dae::HDAE, h::UInt) = hash(dae.d, hash(dae.m, hash(dae.n,
                        hash(dae.v, hash(dae.f, hash(dae.u, hash(dae.g,
                        hash(dae.u̅, hash(dae.g̅, hash(dae.ϕ, hash(dae.ψ,
                        hash(dae.h, hash(dae.v̄, hash(dae.f̄,
                        hash(dae.t₀, hash(dae.q₀, hash(dae.p₀, h)))))))))))))))))

Base.:(==)(dae1::HDAE, dae2::HDAE) = (
                                dae1.d == dae2.d
                             && dae1.m == dae2.m
                             && dae1.n == dae2.n
                             && dae1.v == dae2.v
                             && dae1.f == dae2.f
                             && dae1.u == dae2.u
                             && dae1.g == dae2.g
                             && dae1.u̅ == dae2.u̅
                             && dae1.g̅ == dae2.g̅
                             && dae1.ϕ == dae2.ϕ
                             && dae1.ψ == dae2.ψ
                             && dae1.h == dae2.h
                             && dae1.v̄ == dae2.v̄
                             && dae1.f̄ == dae2.f̄
                             && dae1.t₀ == dae2.t₀
                             && dae1.q₀ == dae2.q₀
                             && dae1.p₀ == dae2.p₀)

function Base.similar(dae::HDAE, q₀, p₀, λ₀=get_λ₀(q₀, dae.λ₀); kwargs...)
    similar(dae, dae.t₀, q₀, p₀, λ₀; kwargs...)
end

function Base.similar(dae::HDAE, t₀::TT, q₀::AbstractArray{DT}, p₀::AbstractArray{DT}, λ₀::AbstractArray{DT}=get_λ₀(q₀, dae.λ₀);
                      v̄=dae.v̄, f̄=dae.f̄, parameters=dae.parameters, periodicity=dae.periodicity) where {DT  <: Number, TT <: Number}
    @assert dae.d == size(q₀,1) == size(p₀,1)
    @assert dae.m == size(λ₀,1)
    HDAE(dae.v, dae.f, dae.u, dae.g, dae.u̅, dae.g̅, dae.ϕ, dae.ψ, dae.h, t₀, q₀, p₀, λ₀;
         v̄=v̄, f̄=f̄, parameters=parameters, periodicity=periodicity)
end

@inline Base.ndims(equation::HDAE) = equation.d
@inline CommonFunctions.nconstraints(equation::HDAE) = equation.m
@inline CommonFunctions.periodicity(equation::HDAE) = equation.periodicity

function get_function_tuple(equation::HDAE{DT,TT,VT,FT,UT,GT,U̅T,G̅T,ϕT,ψT,HT,V̄T,F̄T,Nothing}) where {DT, TT, VT, FT, UT, GT, U̅T, G̅T, ϕT, ψT, HT, V̄T, F̄T}
    NamedTuple{(:v, :f, :u, :g, :u̅, :g̅, :ϕ, :ψ, :h, :v̄, :f̄)}((
        equation.v, equation.f, equation.u, equation.g,
        equation.u̅, equation.g̅, equation.ϕ, equation.ψ,
        equation.h, equation.v̄, equation.f̄))
end

function get_function_tuple(equation::HDAE{DT,TT,VT,FT,UT,GT,U̅T,G̅T,ϕT,ψT,HT,V̄T,F̄T,PT}) where {DT, TT, VT, FT, UT, GT, U̅T, G̅T, ϕT, ψT, HT, V̄T, F̄T, PT <: NamedTuple}
    vₚ = (t,q,p,v)   -> equation.v(t, q, p, v, equation.parameters)
    fₚ = (t,q,p,f)   -> equation.f(t, q, p, f, equation.parameters)
    uₚ = (t,q,p,λ,u) -> equation.u(t, q, p, λ, u, equation.parameters)
    gₚ = (t,q,p,λ,g) -> equation.g(t, q, p, λ, g, equation.parameters)
    u̅ₚ = (t,q,p,λ,u̅) -> equation.u̅(t, q, p, λ, u̅, equation.parameters)
    g̅ₚ = (t,q,p,λ,g̅) -> equation.g̅(t, q, p, λ, g̅, equation.parameters)
    ϕₚ = (t,q,p,ϕ)   -> equation.ϕ(t, q, p, ϕ, equation.parameters)
    ψₚ = (t,q,p,v,f,ψ) -> equation.ψ(t, q, p, v, f, ψ, equation.parameters)
    hₚ = (t,q,p) -> equation.h(t, q, p, equation.parameters)
    v̄ₚ = (t,q,p,v) -> equation.v̄(t, q, p, v, equation.parameters)
    f̄ₚ = (t,q,p,f) -> equation.f̄(t, q, p, f, equation.parameters)

    names = (:v, :f, :u, :g, :u̅, :g̅, :ϕ, :ψ, :h, :v̄, :f̄)
    equs  = (vₚ, fₚ, uₚ, gₚ, u̅ₚ, g̅ₚ, ϕₚ, ψₚ, hₚ, v̄ₚ, f̄ₚ)

    NamedTuple{names}(equs)
end
