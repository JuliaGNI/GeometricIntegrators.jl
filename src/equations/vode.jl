@doc raw"""
`VODE`: Variational Ordinary Differential Equation *EXPERIMENTAL*

Defines an implicit initial value problem
```math
\begin{aligned}
\dot{q} (t) &= v(t) , &
q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), v(t)) , &
p(t_{0}) &= p_{0} , \\
p(t) &= ϑ(t, q(t), v(t))
\end{aligned}
```
with vector field ``f``, the momentum defined by ``p``, initial conditions ``(q_{0}, p_{0})`` and the solution
``(q,p)`` taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}``.
This is a special case of a differential algebraic equation with dynamical
variables ``(q,p)`` and algebraic variable ``v``.

### Fields

* `d`: dimension of dynamical variables ``q`` and ``p`` as well as the vector fields ``f`` and ``p``
* `ϑ`: function determining the momentum
* `f`: function computing the vector field
* `g`: function determining the projection, given by ∇ϑ(q)⋅λ
* `v̄`: function computing an initial guess for the velocity field ``v``` (optional)
* `f̄`: function computing an initial guess for the force field ``f`` (optional)
* `h`: function computing the Hamiltonian (optional)
* `Ω`: symplectic matrix (optional)
* `∇H`: gradient of the Hamiltonian (optional)
* `t₀`: initial time (optional)
* `q₀`: initial condition for `q`
* `p₀`: initial condition for `p`
* `λ₀`: initial condition for `λ` (optional)

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
struct VODE{dType <: Number, tType <: Number, ϑType <: Function, fType <: Function, gType <: Function,
            v̄Type <: Function, f̄Type <: Function, hType <: Union{Function,Nothing},
            ΩType <: Union{Function,Nothing}, ∇HType <: Union{Function,Nothing},
            pType <: Union{NamedTuple,Nothing}, N} <: AbstractEquationPODE{dType, tType}

    d::Int
    m::Int
    n::Int
    ϑ::ϑType
    f::fType
    g::gType
    v̄::v̄Type
    f̄::f̄Type
    h::hType
    Ω::ΩType
    ∇H::∇HType
    t₀::tType
    q₀::Array{dType, N}
    p₀::Array{dType, N}
    λ₀::Array{dType, N}
    parameters::pType
    periodicity::Vector{dType}

    function VODE(DT::DataType, N::Int, d::Int, n::Int, ϑ::ϑType, f::fType, g::gType,
                  t₀::tType, q₀::AbstractArray{dType}, p₀::AbstractArray{dType}, λ₀::AbstractArray{dType};
                  v̄::v̄Type=(t,q,v)->nothing, f̄::f̄Type=f, h::hType=nothing, Ω::ΩType=nothing, ∇H::∇HType=nothing,
                  parameters::pType=nothing, periodicity=zeros(DT,d)) where {
                        dType <: Number, tType <: Number, ϑType <: Function,
                        fType <: Function, gType <: Function,
                        v̄Type <: Function, f̄Type <: Function,
                        hType <: Union{Function,Nothing},
                        ΩType <: Union{Function,Nothing}, ∇HType <: Union{Function,Nothing},
                        pType <: Union{NamedTuple,Nothing}}

        @assert d == size(q₀,1) == size(p₀,1) == size(λ₀,1)
        @assert n == size(q₀,2) == size(p₀,2)
        @assert ndims(q₀) == ndims(p₀) == N ∈ (1,2)

        new{DT, tType, ϑType, fType, gType, v̄Type, f̄Type, hType, ΩType, ∇HType, pType, N}(d, d, n, ϑ, f, g, v̄, f̄, h, Ω, ∇H,
                t₀, convert(Array{DT}, q₀), convert(Array{DT}, p₀), convert(Array{DT}, λ₀), parameters, periodicity)
    end
end

function VODE(ϑ, f, g, t₀::Number, q₀::AbstractArray{DT}, p₀::AbstractArray{DT}, λ₀::AbstractArray{DT}=zero(q₀); kwargs...) where {DT}
    VODE(DT, ndims(q₀), size(q₀,1), size(q₀,2), ϑ, f, g, t₀, q₀, p₀, λ₀; kwargs...)
end

function VODE(ϑ, f, g, q₀::AbstractArray, p₀::AbstractArray, λ₀::AbstractArray=zero(q₀); kwargs...)
    VODE(ϑ, f, g, zero(eltype(q₀)), q₀, p₀, λ₀; kwargs...)
end

Base.hash(ode::VODE, h::UInt) = hash(ode.d, hash(ode.m, hash(ode.n,
          hash(ode.ϑ, hash(ode.f, hash(ode.g, hash(ode.v̄, hash(ode.f̄,
          hash(ode.h, hash(ode.Ω, hash(ode.∇H,
          hash(ode.t₀, hash(ode.q₀, hash(ode.p₀, hash(ode.λ₀,
          hash(ode.parameters, hash(ode.periodicity, h)))))))))))))))))

Base.:(==)(ode1::VODE, ode2::VODE) = (
                                ode1.d == ode2.d
                             && ode1.m == ode2.m
                             && ode1.n == ode2.n
                             && ode1.ϑ == ode2.ϑ
                             && ode1.f == ode2.f
                             && ode1.g == ode2.g
                             && ode1.v̄ == ode2.v̄
                             && ode1.f̄ == ode2.f̄
                             && ode1.h == ode2.h
                             && ode1.Ω == ode2.Ω
                             && ode1.∇H == ode2.∇H
                             && ode1.t₀ == ode2.t₀
                             && ode1.q₀ == ode2.q₀
                             && ode1.p₀ == ode2.p₀
                             && ode1.λ₀ == ode2.λ₀
                             && ode1.parameters == ode2.parameters
                             && ode1.periodicity == ode2.periodicity)

function Base.similar(ode::VODE, q₀, p₀, λ₀=get_λ₀(q₀, ode.λ₀); kwargs...)
    similar(ode, ode.t₀, q₀, p₀, λ₀; kwargs...)
end

function Base.similar(ode::VODE, t₀::TT, q₀::AbstractArray{DT}, p₀::AbstractArray{DT}, λ₀::AbstractArray{DT}=get_λ₀(q₀, ode.λ₀);
                      v̄=ode.v̄, f̄=ode.f̄, h=ode.h, Ω=ode.Ω, ∇H=ode.∇H, parameters=ode.parameters, periodicity=ode.periodicity) where {DT  <: Number, TT <: Number}
    @assert ode.d == size(q₀,1) == size(p₀,1) == size(λ₀,1)
    VODE(ode.ϑ, ode.f, ode.g, t₀, q₀, p₀, λ₀; v̄=v̄, f̄=f̄, h=h, Ω=Ω, ∇H=∇H,
         parameters=parameters, periodicity=periodicity)
end

@inline Base.ndims(ode::VODE) = ode.d

@inline CommonFunctions.periodicity(equation::VODE) = equation.periodicity

function get_function_tuple(equation::VODE{DT,TT,ϑT,FT,GT,V̄T,F̄T,HT,ΩT,∇HT,Nothing}) where {DT, TT, ϑT, FT, GT, V̄T, F̄T, HT, ΩT, ∇HT}
    names = (:ϑ,:f,:g,:v̄,:f̄)
    equs  = (equation.ϑ, equation.f, equation.g, equation.v̄, equation.f̄)

    if HT != Nothing
        names = (names..., :h)
        equs  = (equs..., equation.h)
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

function get_function_tuple(equation::VODE{DT,TT,ϑT,FT,GT,V̄T,F̄T,HT,ΩT,∇HT,PT}) where {DT, TT, ϑT, FT, GT, V̄T, F̄T, HT, ΩT, ∇HT, PT <: NamedTuple}
    ϑₚ = (t,q,v,ϑ) -> equation.ϑ(t, q, v, ϑ, equation.parameters)
    fₚ = (t,q,v,f) -> equation.f(t, q, v, f, equation.parameters)
    gₚ = (t,q,v,g) -> equation.g(t, q, v, g, equation.parameters)
    v̄ₚ = (t,q,v)   -> equation.v̄(t, q, v, equation.parameters)
    f̄ₚ = (t,q,v,f) -> equation.f̄(t, q, v, f, equation.parameters)

    names = (:ϑ, :f, :g, :v̄, :f̄)
    equs  = (ϑₚ, fₚ, gₚ, v̄ₚ, f̄ₚ)

    if HT != Nothing
        hₚ = (t,q) -> equation.h(t, q, equation.parameters)
        names = (names..., :h)
        equs  = (equs..., hₚ)
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
