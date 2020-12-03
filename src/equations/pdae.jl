@doc raw"""
`PDAE`: Partitioned Differential Algebraic Equation

Defines a partitioned differential algebraic initial value problem
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t), p(t)) + u(t, q(t), p(t), \lambda(t)) , & q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), p(t)) + r(t, q(t), p(t), \lambda(t)) , & p(t_{0}) &= p_{0} , \\
0 &= \phi (t, q(t), p(t), \lambda(t)) , & \lambda(t_{0}) &= \lambda_{0} ,
\end{aligned}
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
* `u`: function computing the projection for ``q``
* `g`: function computing the projection for ``p``
* `ϕ`: algebraic constraint
* `v̄`: function computing an initial guess for the velocity field ``v`` (optional)
* `f̄`: function computing an initial guess for the force field ``f`` (optional)
* `h`: function computing the Hamiltonian (optional)
* `t₀`: initial time
* `q₀`: initial condition for dynamical variable ``q``
* `p₀`: initial condition for dynamical variable ``p``
* `λ₀`: initial condition for algebraic variable ``\lambda``

"""
struct PDAE{dType <: Number, tType <: Number,
            vType <: Function, fType <: Function,
            uType <: Function, gType <: Function,
            ϕType <: Function,
            v̄Type <: Function, f̄Type <: Function,
            hType <: Union{Function,Nothing},
            pType <: Union{NamedTuple,Nothing}, N} <: AbstractEquationPDAE{dType, tType}

    d::Int
    m::Int
    n::Int
    v::vType
    f::fType
    u::uType
    g::gType
    ϕ::ϕType
    v̄::v̄Type
    f̄::f̄Type
    h::hType
    t₀::tType
    q₀::Array{dType, N}
    p₀::Array{dType, N}
    λ₀::Array{dType, N}
    parameters::pType
    periodicity::Vector{dType}

    function PDAE(DT::DataType, N::Int, d::Int, m::Int, n::Int,
                  v::vType, f::fType, u::uType, g::gType, ϕ::ϕType,
                  t₀::tType, q₀::AbstractArray{dType}, p₀::AbstractArray{dType}, λ₀::AbstractArray{dType};
                  v̄::v̄Type=v, f̄::f̄Type=f, h::hType=nothing, parameters::pType=nothing, periodicity=zeros(DT,d)) where {
                         dType <: Number, tType <: Number,
                         vType <: Function, fType <: Function,
                         uType <: Function, gType <: Function,
                         ϕType <: Function,
                         v̄Type <: Function, f̄Type <: Function,
                         hType <: Union{Function,Nothing},
                         pType <: Union{NamedTuple,Nothing}}

        @assert d == size(q₀,1) == size(p₀,1)
        @assert m == size(λ₀,1)
        @assert n == size(q₀,2) == size(p₀,2) == size(λ₀,2)
        @assert 2n ≥ m
        @assert ndims(q₀) == ndims(p₀) == ndims(λ₀) == N ∈ (1,2)

        new{DT, tType, vType, fType, uType, gType, ϕType, v̄Type, f̄Type, hType, pType, N}(d, m, n, v, f, u, g, ϕ, v̄, f̄, h, t₀,
                convert(Array{DT}, q₀), convert(Array{DT}, p₀), convert(Array{DT}, λ₀),
                parameters, periodicity)
    end
end

function PDAE(v, f, u, g, ϕ, t₀, q₀::AbstractArray{DT}, p₀::AbstractArray{DT}, λ₀::AbstractArray{DT}; kwargs...) where {DT}
    PDAE(DT, ndims(q₀), size(q₀,1), size(λ₀,1), size(q₀,2), v, f, u, g, ϕ, t₀, q₀, p₀, λ₀; kwargs...)
end

function PDAE(v, f, u, g, ϕ, q₀, p₀, λ₀; kwargs...)
    PDAE(v, f, u, g, ϕ, zero(eltype(q₀)), q₀, p₀, λ₀; kwargs...)
end

Base.hash(dae::PDAE, h::UInt) = hash(dae.d, hash(dae.m, hash(dae.n,
        hash(dae.v, hash(dae.f, hash(dae.u, hash(dae.g, 
        hash(dae.v̄, hash(dae.f̄, hash(dae.h,
        hash(dae.t₀, hash(dae.q₀, hash(dae.p₀, hash(dae.λ₀,
        hash(dae.periodicity, hash(dae.parameters, h))))))))))))))))

Base.:(==)(dae1::PDAE, dae2::PDAE) = (
                                dae1.d == dae2.d
                             && dae1.m == dae2.m
                             && dae1.n == dae2.n
                             && dae1.v == dae2.v
                             && dae1.f == dae2.f
                             && dae1.u == dae2.u
                             && dae1.g == dae2.g
                             && dae1.ϕ == dae2.ϕ
                             && dae1.v̄ == dae2.v̄
                             && dae1.f̄ == dae2.f̄
                             && dae1.h == dae2.h
                             && dae1.t₀ == dae2.t₀
                             && dae1.q₀ == dae2.q₀
                             && dae1.p₀ == dae2.p₀
                             && dae1.λ₀ == dae2.λ₀
                             && dae1.parameters == dae1.parameters
                             && dae1.periodicity == dae1.periodicity)

function Base.similar(dae::PDAE, q₀, p₀, λ₀=get_λ₀(q₀, dae.λ₀); kwargs...)
    similar(dae, dae.t₀, q₀, p₀, λ₀; kwargs...)
end

function Base.similar(dae::PDAE, t₀::TT, q₀::AbstractArray{DT}, p₀::AbstractArray{DT}, λ₀::AbstractArray{DT}=get_λ₀(q₀, dae.λ₀);
                      v̄=dae.v̄, f̄=dae.f̄, h=dae.h, parameters=dae.parameters, periodicity=dae.periodicity) where {DT  <: Number, TT <: Number}
    @assert dae.d == size(q₀,1) == size(p₀,1)
    @assert dae.m == size(λ₀,1)
    PDAE(dae.v, dae.f, dae.u, dae.g, dae.ϕ, t₀, q₀, p₀, λ₀; v̄=v̄, f̄=f̄, h=h, parameters=parameters, periodicity=periodicity)
end

@inline Base.ndims(equation::PDAE) = equation.d
@inline CommonFunctions.nconstraints(equation::PDAE) = equation.m
@inline CommonFunctions.periodicity(equation::PDAE) = equation.periodicity

function get_function_tuple(equation::PDAE{DT,TT,VT,FT,UT,GT,ϕT,V̄T,F̄T,HT,Nothing}) where {DT, TT, VT, FT, UT, GT, ϕT, V̄T, F̄T, HT}
    names = (:v,:f,:u,:g,:ϕ,:v̄,:f̄)
    equs  = (equation.v, equation.f, equation.u, equation.g, equation.ϕ, equation.v̄, equation.f̄)

    if HT != Nothing
        names = (names..., :h)
        equs  = (equs..., equation.h)
    end

    NamedTuple{names}(equs)
end

function get_function_tuple(equation::PDAE{DT,TT,VT,FT,UT,GT,ϕT,V̄T,F̄T,HT,PT}) where {DT, TT, VT, FT, UT, GT, ϕT, V̄T, F̄T, HT, PT <: NamedTuple}
    vₚ = (t,q,p,v) -> equation.v(t, q, p, v, equation.parameters)
    fₚ = (t,q,p,f) -> equation.f(t, q, p, f, equation.parameters)
    uₚ = (t,q,p,λ,u) -> equation.u(t, q, p, λ, u, equation.parameters)
    gₚ = (t,q,p,λ,g) -> equation.g(t, q, p, λ, g, equation.parameters)
    ϕₚ = (t,q,p,ϕ) -> equation.ϕ(t, q, p, ϕ, equation.parameters)
    v̄ₚ = (t,q,p,v) -> equation.v̄(t, q, p, v, equation.parameters)
    f̄ₚ = (t,q,p,f) -> equation.f̄(t, q, p, f, equation.parameters)

    names = (:v, :f, :u, :g, :ϕ, :v̄, :f̄)
    equs  = (vₚ, fₚ, uₚ, gₚ, ϕₚ, v̄ₚ, f̄ₚ)

    if HT != Nothing
        hₚ = (t,q,p) -> equation.h(t, q, p, equation.parameters)
        names = (names..., :h)
        equs  = (equs..., hₚ)
    end

    NamedTuple{names}(equs)
end
