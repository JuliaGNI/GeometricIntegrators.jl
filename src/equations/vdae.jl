@doc raw"""
`VDAE`: Variational Differential Algebraic Equation *EXPERIMENTAL*

Defines an implicit initial value problem
```math
\begin{align*}
\dot{q} (t) &= v(t) + λ(t), &
q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), v(t)) + g(t, q(t), λ(t)) + \obar{g}(t, q(t), μ(t)) , &
p(t_{0}) &= p_{0} , \\
p(t) &= ϑ(t, q(t), v(t)) , \\
0 &= \phi (t, q(t), p(t)) , \\
0 &= \psi (t, q(t), p(t), \dot{q}(t), \dot{p}(t)) ,
\end{align*}
```
with vector field ``f``, the momentum defined by ``p``, initial conditions ``(q_{0}, p_{0})`` and the solution
``(q,p)`` taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}`` and
the algebraic variables ``(v, \lambda, \mu)`` taking values in
``\mathbb{R}^{d} \times \mathbb{R}^{d} \times \mathbb{R}^{m}``.
This is a special case of a differential algebraic equation with dynamical
variables ``(q,p)`` and algebraic variables ``v``, ``\lambda`` and ``\mu``.

### Fields

* `d`: dimension of dynamical variables ``q`` and ``p`` as well as the vector fields ``f`` and ``p``
* `ϑ`: function determining the momentum
* `f`: function computing the vector field
* `g`: function determining the primary projection, usually given by ∇ϑ(q)⋅λ
* `g̅`: function determining the secondary projection, usually given by λ⋅∇ϑ(q)
* `ϕ`: primary constraints, usually given by p-ϑ(q)
* `ψ`: secondary constraints, usually given by ṗ-q̇⋅∇ϑ(q)
* `v`: function computing an initial guess for the velocity field (optional)
* `Ω`: symplectic matrix (optional)
* `∇H`: gradient of the Hamiltonian (optional)
* `t₀`: initial time (optional)
* `q₀`: initial condition for `q`
* `p₀`: initial condition for `p`
* `λ₀`: initial condition for `λ` (optional)
* `μ₀`: initial condition for `μ` (optional)

The functions `ϑ` and `f` must have the interface
```julia
    function ϑ(t, q, v, p)
        p[1] = ...
        p[2] = ...
        ...
    end
```
and
```julia
    function f(t, q, v, f)
        f[1] = ...
        f[2] = ...
        ...
    end
```
where `t` is the current time, `q` is the current solution vector, `v` is the
current velocity and `f` and `p` are the vectors which hold the result of
evaluating the functions ``f`` and ``ϑ`` on `t`, `q` and `v`.
The funtions `g` and `v` are specified by
```julia
    function g(t, q, λ, g)
        g[1] = ...
        g[2] = ...
        ...
    end
```
and
```julia
    function v(t, q, p, v)
        v[1] = ...
        v[2] = ...
        ...
    end
```
"""
struct VDAE{dType <: Number, tType <: Number, ϑType <: Function,
            fType <: Function, gType <: Function, g̅Type <: Function,
            ϕType <: Function, ψType <: Function, vType <: Function,
            ΩType <: Union{Function,Nothing}, ∇HType <: Union{Function,Nothing},
            pType <: Union{Tuple,Nothing}, N} <: Equation{dType, tType}

    d::Int
    m::Int
    n::Int
    ϑ::ϑType
    f::fType
    g::gType
    g̅::g̅Type
    ϕ::ϕType
    ψ::ψType
    v::vType
    Ω::ΩType
    ∇H::∇HType
    t₀::tType
    q₀::Array{dType, N}
    p₀::Array{dType, N}
    λ₀::Array{dType, N}
    μ₀::Array{dType, N}
    parameters::pType
    periodicity::Vector{dType}

    function VDAE(DT::DataType, N::Int, d::Int, m::Int, n::Int, ϑ::ϑType, f::fType,
                  g::gType, g̅::g̅Type, ϕ::ϕType, ψ::ψType, t₀::tType,
                  q₀::DenseArray{dType}, p₀::DenseArray{dType},
                  λ₀::DenseArray{dType}, μ₀::DenseArray{dType};
                  v::vType=nothing, Ω::ΩType=nothing, ∇H::∇HType=nothing,
                  parameters=nothing, periodicity=zeros(DT,d)) where {
                        dType <: Number, tType <: Number, ϑType <: Function,
                        fType <: Function, gType <: Function, g̅Type <: Function,
                        ϕType <: Function, ψType <: Function, vType <: Function,
                        ΩType <: Union{Function,Nothing}, ∇HType <: Union{Function,Nothing}}

        @assert d == size(q₀,1) == size(p₀,1) == size(λ₀,1)
        @assert m == size(μ₀,1)
        @assert n == size(q₀,2) == size(p₀,2)
        @assert ndims(q₀) == ndims(p₀) == N ∈ (1,2)

        new{DT, tType, ϑType, fType, gType, g̅Type, ϕType, ψType, vType, ΩType, ∇HType, typeof(parameters), N}(d, m, n, ϑ, f, g, g̅, ϕ, ψ, v, Ω, ∇H, t₀,
                convert(Array{DT}, q₀), convert(Array{DT}, p₀), convert(Array{DT}, λ₀), convert(Array{DT}, μ₀), parameters, periodicity)
    end
end


function VDAE(ϑ, f, g, g̅, ϕ, ψ, t₀::Number, q₀::DenseArray{DT}, p₀::DenseArray{DT}, λ₀::DenseArray{DT}=zero(q₀), μ₀::DenseArray{DT}=zero(q₀); kwargs...) where {DT}
    VDAE(DT, ndims(q₀), size(q₀,1), size(μ₀,1), size(q₀,2), ϑ, f, g, g̅, ϕ, ψ, t₀, q₀, p₀, λ₀, μ₀; kwargs...)
end

function VDAE(ϑ, f, g, g̅, ϕ, ψ, q₀::DenseArray, p₀::DenseArray, λ₀::DenseArray=zero(q₀), μ₀::DenseArray=zero(q₀); kwargs...)
    VDAE(ϑ, f, g, g̅, ϕ, ψ, zero(eltype(q₀)), q₀, p₀, λ₀, μ₀; kwargs...)
end

Base.hash(ode::VDAE, h::UInt) = hash(ode.d, hash(ode.m, hash(ode.n,
          hash(ode.ϑ, hash(ode.f, hash(ode.g, hash(ode.g̅, hash(ode.ϕ,
          hash(ode.ψ, hash(ode.v, hash(ode.Ω, hash(ode.∇H,
          hash(ode.t₀, hash(ode.q₀, hash(ode.p₀, hash(ode.λ₀, hash(ode.μ₀,
          hash(ode.parameters, hash(ode.periodicity, h)))))))))))))))))))

Base.:(==)(ode1::VDAE, ode2::VDAE) = (
                                ode1.d == ode2.d
                             && ode1.m == ode2.m
                             && ode1.n == ode2.n
                             && ode1.ϑ == ode2.ϑ
                             && ode1.f == ode2.f
                             && ode1.g == ode2.g
                             && ode1.g̅ == ode2.g̅
                             && ode1.ϕ == ode2.ϕ
                             && ode1.ψ == ode2.ψ
                             && ode1.v == ode2.v
                             && ode1.Ω == ode2.Ω
                             && ode1.∇H == ode2.∇H
                             && ode1.t₀ == ode2.t₀
                             && ode1.q₀ == ode2.q₀
                             && ode1.p₀ == ode2.p₀
                             && ode1.λ₀ == ode2.λ₀
                             && ode1.μ₀ == ode2.μ₀
                             && ode1.parameters == ode2.parameters
                             && ode1.periodicity == ode2.periodicity)

function Base.similar(ode::VDAE, q₀, p₀, λ₀=get_λ₀(q₀, ode.λ₀), μ₀=get_λ₀(q₀, ode.μ₀); kwargs...)
    similar(ode, ode.t₀, q₀, p₀, λ₀, μ₀; kwargs...)
end

function Base.similar(ode::VDAE, t₀::TT, q₀::DenseArray{DT}, p₀::DenseArray{DT},
                      λ₀::DenseArray{DT}=get_λ₀(q₀, ode.λ₀), μ₀=get_λ₀(q₀, ode.μ₀);
                      parameters=ode.parameters, periodicity=ode.periodicity) where {DT  <: Number, TT <: Number}
    @assert ode.d == size(q₀,1) == size(p₀,1) == size(λ₀,1)
    VDAE(ode.ϑ, ode.f, ode.g, ode.g̅, ode.ϕ, ode.ψ, t₀, q₀, p₀, λ₀, μ₀;
         v=ode.v, Ω=ode.Ω, ∇H=ode.∇H, parameters=parameters, periodicity=periodicity)
end

Base.ndims(ode::VDAE) = ode.d
