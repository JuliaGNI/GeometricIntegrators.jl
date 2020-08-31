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
* `h`: function computing the Hamiltonian (optional)
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
            ϕType <: Function, ψType <: Function,
            hType <: Union{Function,Nothing}, vType <: Union{Function,Nothing},
            ΩType <: Union{Function,Nothing}, ∇HType <: Union{Function,Nothing},
            pType <: Union{NamedTuple,Nothing}, N} <: AbstractEquationPDAE{dType, tType}

    d::Int
    m::Int
    n::Int
    ϑ::ϑType
    f::fType
    g::gType
    g̅::g̅Type
    ϕ::ϕType
    ψ::ψType
    h::hType
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
                  q₀::AbstractArray{dType}, p₀::AbstractArray{dType},
                  λ₀::AbstractArray{dType}, μ₀::AbstractArray{dType};
                  h::hType=nothing, v::vType=nothing, Ω::ΩType=nothing, ∇H::∇HType=nothing,
                  parameters::pType=nothing, periodicity=zeros(DT,d)) where {
                        dType <: Number, tType <: Number,
                        ϑType <: Function, fType <: Function, gType <: Function,
                        g̅Type <: Function, ϕType <: Function, ψType <: Function,
                        hType <: Union{Function,Nothing}, vType <: Union{Function,Nothing},
                        ΩType <: Union{Function,Nothing}, ∇HType <: Union{Function,Nothing},
                        pType <: Union{NamedTuple,Nothing}}

        @assert d == size(q₀,1) == size(p₀,1) == size(λ₀,1)
        @assert m == size(μ₀,1)
        @assert n == size(q₀,2) == size(p₀,2)
        @assert ndims(q₀) == ndims(p₀) == N ∈ (1,2)

        new{DT, tType, ϑType, fType, gType, g̅Type, ϕType, ψType, hType, vType, ΩType, ∇HType, pType, N}(d, m, n, ϑ, f, g, g̅, ϕ, ψ, h, v, Ω, ∇H, t₀,
                convert(Array{DT}, q₀), convert(Array{DT}, p₀), convert(Array{DT}, λ₀), convert(Array{DT}, μ₀), parameters, periodicity)
    end
end


function VDAE(ϑ, f, g, g̅, ϕ, ψ, t₀::Number, q₀::AbstractArray{DT}, p₀::AbstractArray{DT}, λ₀::AbstractArray{DT}=zero(q₀), μ₀::AbstractArray{DT}=zero(q₀); kwargs...) where {DT}
    VDAE(DT, ndims(q₀), size(q₀,1), size(μ₀,1), size(q₀,2), ϑ, f, g, g̅, ϕ, ψ, t₀, q₀, p₀, λ₀, μ₀; kwargs...)
end

function VDAE(ϑ, f, g, g̅, ϕ, ψ, q₀::AbstractArray, p₀::AbstractArray, λ₀::AbstractArray=zero(q₀), μ₀::AbstractArray=zero(q₀); kwargs...)
    VDAE(ϑ, f, g, g̅, ϕ, ψ, zero(eltype(q₀)), q₀, p₀, λ₀, μ₀; kwargs...)
end

Base.hash(dae::VDAE, h::UInt) = hash(dae.d, hash(dae.m, hash(dae.n,
          hash(dae.ϑ, hash(dae.f, hash(dae.g, hash(dae.g̅, hash(dae.ϕ,
          hash(dae.ψ, hash(dae.h, hash(dae.v, hash(dae.Ω, hash(dae.∇H,
          hash(dae.t₀, hash(dae.q₀, hash(dae.p₀, hash(dae.λ₀, hash(dae.μ₀,
          hash(dae.parameters, hash(dae.periodicity, h))))))))))))))))))))

Base.:(==)(dae1::VDAE, dae2::VDAE) = (
                                dae1.d == dae2.d
                             && dae1.m == dae2.m
                             && dae1.n == dae2.n
                             && dae1.ϑ == dae2.ϑ
                             && dae1.f == dae2.f
                             && dae1.g == dae2.g
                             && dae1.g̅ == dae2.g̅
                             && dae1.ϕ == dae2.ϕ
                             && dae1.ψ == dae2.ψ
                             && dae1.h == dae2.h
                             && dae1.v == dae2.v
                             && dae1.Ω == dae2.Ω
                             && dae1.∇H == dae2.∇H
                             && dae1.t₀ == dae2.t₀
                             && dae1.q₀ == dae2.q₀
                             && dae1.p₀ == dae2.p₀
                             && dae1.λ₀ == dae2.λ₀
                             && dae1.μ₀ == dae2.μ₀
                             && dae1.parameters == dae2.parameters
                             && dae1.periodicity == dae2.periodicity)

function Base.similar(dae::VDAE, q₀, p₀, λ₀=get_λ₀(q₀, dae.λ₀), μ₀=get_λ₀(q₀, dae.μ₀); kwargs...)
    similar(dae, dae.t₀, q₀, p₀, λ₀, μ₀; kwargs...)
end

function Base.similar(dae::VDAE, t₀::TT, q₀::AbstractArray{DT}, p₀::AbstractArray{DT},
                      λ₀::AbstractArray{DT}=get_λ₀(q₀, dae.λ₀), μ₀=get_λ₀(q₀, dae.μ₀);
                      h=dae.h, v=dae.v, Ω=dae.Ω, ∇H=dae.∇H, parameters=dae.parameters, periodicity=dae.periodicity) where {DT  <: Number, TT <: Number}
    @assert dae.d == size(q₀,1) == size(p₀,1) == size(λ₀,1)
    VDAE(dae.ϑ, dae.f, dae.g, dae.g̅, dae.ϕ, dae.ψ, t₀, q₀, p₀, λ₀, μ₀;
         h=h, v=v, Ω=Ω, ∇H=∇H, parameters=parameters, periodicity=periodicity)
end

@inline Base.ndims(equation::VDAE) = equation.d
@inline CommonFunctions.nconstraints(equation::VDAE) = equation.m
@inline CommonFunctions.periodicity(equation::VDAE) = equation.periodicity


function get_function_tuple(equation::VDAE{DT,TT,θT,FT,GT,G̅T,ϕT,ψT,HT,VT,ΩT,∇HT,Nothing}) where {DT,TT,θT,FT,GT,G̅T,ϕT,ψT,HT,VT,ΩT,∇HT}
    names = (:ϑ,:f,:g,:g̅,:ϕ,:ψ)
    equs  = (equation.ϑ, equation.f, equation.g, equation.g̅, equation.ϕ, equation.ψ)

    if HT != Nothing
        names = (names..., :h)
        equs  = (equs..., equation.h)
    end

    if VT != Nothing
        names = (names..., :v)
        equs  = (equs..., equation.v)
    end

    if ΩT != Nothing
        names = (names..., :Ω)
        equs  = (equs..., equation.Ω)
    end

    if ∇HT != Nothing
        names = (names..., :∇H)
        equs  = (equs..., equation.∇H)
    end

    NamedTuple{names}(equs)
end

function get_function_tuple(equation::VDAE{DT,TT,θT,FT,GT,G̅T,ϕT,ψT,HT,VT,ΩT,∇HT,PT}) where {DT, TT, θT, FT, GT, G̅T, ϕT, ψT, HT, VT, ΩT, ∇HT, PT <: NamedTuple}
    ϑₚ = (t,q,v,ϑ) -> equation.ϑ(t, q, v, ϑ, equation.parameters)
    fₚ = (t,q,v,f) -> equation.f(t, q, v, f, equation.parameters)
    gₚ = (t,q,λ,g) -> equation.g(t, q, λ, g, equation.parameters)
    g̅ₚ = (t,q,λ,g̅) -> equation.g̅(t, q, λ, g̅, equation.parameters)
    ϕₚ = (t,q,v,ϕ) -> equation.ϕ(t, q, v, ϕ, equation.parameters)
    ψₚ = (t,q,v,p,f,ψ) -> equation.ψ(t, q, v, p, f, ψ, equation.parameters)

    names = (:ϑ, :f, :g, :g̅, :ϕ, :ψ)
    equs  = (ϑₚ, fₚ, gₚ, g̅ₚ, ϕₚ, ψₚ)

    if HT != Nothing
        hₚ = (t,q) -> equation.h(t, q, equation.parameters)
        names = (names..., :h)
        equs  = (equs..., hₚ)
    end

    if VT != Nothing
        vₚ = (t,q,v) -> equation.v(t, q, v, equation.parameters)
        names = (names..., :v)
        equs  = (equs..., vₚ)
    end

    if ΩT != Nothing
        Ωₚ = (t,q,Ω) -> equation.Ω(t, q, Ω, equation.parameters)
        names = (names..., :Ω)
        equs  = (equs..., Ωₚ)
    end

    if ∇HT != Nothing
        ∇Hₚ = (t,q,∇H) -> equation.∇H(t, q, ∇H, equation.parameters)
        names = (names..., :∇H)
        equs  = (equs..., ∇Hₚ)
    end

    NamedTuple{names}(equs)
end
