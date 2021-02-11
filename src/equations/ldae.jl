@doc raw"""
`LDAE`: Variational Differential Algebraic Equation *EXPERIMENTAL*

Defines an implicit initial value problem
```math
\begin{aligned}
\dot{q} (t) &= v(t) + \lambda(t), &
q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), v(t)) + g(t, q(t), \lambda(t)) + \bar{g} (t, q(t), \mu(t)) , &
p(t_{0}) &= p_{0} , \\
p(t) &= ϑ(t, q(t), v(t)) , \\
0 &= \phi (t, q(t), p(t)) , \\
0 &= \psi (t, q(t), p(t), \dot{q}(t), \dot{p}(t)) ,
\end{aligned}
```
with vector field ``f``, the momentum defined by ``p``, initial conditions ``(q_{0}, p_{0})`` and the solution
``(q,p)`` taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}`` and
the algebraic variables ``(v, \lambda, \mu)`` taking values in
``\mathbb{R}^{d} \times \mathbb{R}^{d} \times \mathbb{R}^{m}``.
This is a special case of a differential algebraic equation with dynamical
variables ``(q,p)`` and algebraic variables ``v``, ``\lambda`` and ``\mu``.

### Parameters

* `DT <: Number`: data type
* `TT <: Real`: time step type
* `AT <: AbstractArray{DT}`: array type
* `ϑType <: Function`: type of `ϑ`
* `fType <: Function`: type of `f`
* `uType <: Function`: type of `u`
* `gType <: Function`: type of `g`
* `ϕType <: Function`: type of `ϕ`
* `ψType <: Function`: type of `ψ`
* `v̄Type <: Function`: type of `v̄`
* `f̄Type <: Function`: type of `f̄`
* `hType <: OptionalFunction`: type of `h`
* `ΩType <: OptionalFunction`: type of `Ω`
* `∇HType <: OptionalFunction`: type of `∇H`
* `pType <: Union{NamedTuple,Nothing}`: parameters type

### Fields

* `d`: dimension of dynamical variables ``q`` and ``p`` as well as the vector fields ``f`` and ``p``
* `ϑ`: function determining the momentum
* `f`: function computing the vector field
* `g`: function determining the primary projection, usually given by ``\nabla \vartheta (q) \cdot \lambda``
* `ḡ`: function determining the secondary projection, usually given by ``\lambda \cdot \nabla \vartheta (q)``
* `ϕ`: primary constraints, usually given by ``p - \vartheta (q)``
* `ψ`: secondary constraints, usually given by ``\dot{p} - \dot{q} \cdot \nabla \vartheta (q)``
* `v̄`: function computing an initial guess for the velocity field ``v`` (optional)
* `f̄`: function computing an initial guess for the force field ``f`` (optional)
* `h`: function computing the Hamiltonian (optional)
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
The funtions `g`, `v̄` and `f̄` are specified by
```julia
    function g(t, q, λ, g)
        g[1] = ...
        g[2] = ...
        ...
    end

    function v̄(t, q, v)
        v[1] = ...
        v[2] = ...
        ...
    end

    function f̄(t, q, v, f)
        f[1] = ...
        f[2] = ...
        ...
    end
```

### Constructors

```julia
LDAE(ϑ, f, g, ḡ, ϕ, ψ, t₀, q₀, p₀, λ₀, μ₀; v̄=(t,q,v)->nothing, f̄=f, h=nothing, Ω=nothing, ∇H=nothing, parameters=nothing, periodicity=zero(q₀[begin]))
LDAE(ϑ, f, g, ḡ, ϕ, ψ, q₀::StateVector, p₀::StateVector, λ₀::StateVector=zero(q₀), μ₀::StateVector=zero(q₀); kwargs...) = LDAE(ϑ, f, g, ḡ, ϕ, ψ, 0.0, q₀, p₀, λ₀, μ₀; kwargs...)
LDAE(ϑ, f, g, ḡ, ϕ, ψ, t₀, q₀::State, p₀::State, λ₀::State=zero(q₀), μ₀::State=zero(q₀); kwargs...) = LDAE(ϑ, f, g, ḡ, ϕ, ψ, t₀, [q₀], [p₀], [λ₀], [μ₀]; kwargs...)
LDAE(ϑ, f, g, ḡ, ϕ, ψ, q₀::State, p₀::State, λ₀::State=zero(q₀), μ₀::State=zero(q₀); kwargs...) = LDAE(ϑ, f, g, ḡ, ϕ, ψ, 0.0, q₀, p₀, λ₀, μ₀; kwargs...)
```

"""
struct LDAE{dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
            ϑType <: Function, fType <: Function,
            gType <: Function, ḡType <: Function,
            ϕType <: Function, ψType <: Function,
            v̄Type <: Function, f̄Type <: Function,
            hType <: OptionalFunction,
            ΩType <: OptionalFunction, ∇HType <: OptionalFunction,
            pType <: Union{NamedTuple,Nothing}} <: AbstractEquationPDAE{dType, tType}

    d::Int
    m::Int
    ϑ::ϑType
    f::fType
    g::gType
    ḡ::ḡType
    ϕ::ϕType
    ψ::ψType
    v̄::v̄Type
    f̄::f̄Type
    h::hType
    Ω::ΩType
    ∇H::∇HType
    t₀::tType
    q₀::Vector{arrayType}
    p₀::Vector{arrayType}
    λ₀::Vector{arrayType}
    μ₀::Vector{arrayType}
    parameters::pType
    periodicity::arrayType

    function LDAE(ϑ::ϑType, f::fType, g::gType, ḡ::ḡType, ϕ::ϕType, ψ::ψType, t₀::tType,
                  q₀::Vector{arrayType}, p₀::Vector{arrayType},
                  λ₀::Vector{arrayType}, μ₀::Vector{arrayType};
                  v̄::v̄Type=(t,q,v)->nothing, f̄::f̄Type=f, h::hType=nothing, Ω::ΩType=nothing, ∇H::∇HType=nothing,
                  parameters::pType=nothing, periodicity=zero(q₀[begin])) where {
                        dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
                        ϑType <: Function, fType <: Function,
                        gType <: Function, ḡType <: Function,
                        ϕType <: Function, ψType <: Function,
                        v̄Type <: Function, f̄Type <: Function,
                        hType <: OptionalFunction,
                        ΩType <: OptionalFunction,
                        ∇HType <: OptionalFunction,
                        pType <: Union{NamedTuple,Nothing}}

        d = length(q₀[begin])
        m = length(μ₀[begin])

        @assert 2d ≥ m

        @assert length(q₀) == length(p₀) == length(λ₀) == length(μ₀)

        @assert all(length(q) == d for q in q₀)
        @assert all(length(p) == d for p in p₀)
        @assert all(length(λ) == d for λ in λ₀)
        @assert all(length(μ) == m for μ in μ₀)

        @assert all([ndims(q) == ndims(p) == ndims(λ) == ndims(μ) for (q,p,λ,μ) in zip(q₀,p₀,λ₀,μ₀)])
        
        new{dType, tType, arrayType, ϑType, fType, gType, ḡType, ϕType, ψType, v̄Type, f̄Type, hType, ΩType, ∇HType, pType}(
                d, m, ϑ, f, g, ḡ, ϕ, ψ, v̄, f̄, h, Ω, ∇H, t₀, q₀, p₀, λ₀, μ₀, parameters, periodicity)
    end
end

LDAE(ϑ, f, g, ḡ, ϕ, ψ, q₀::StateVector, p₀::StateVector, λ₀::StateVector=zero(q₀), μ₀::StateVector=zero(q₀); kwargs...) = LDAE(ϑ, f, g, ḡ, ϕ, ψ, 0.0, q₀, p₀, λ₀, μ₀; kwargs...)
LDAE(ϑ, f, g, ḡ, ϕ, ψ, t₀, q₀::State, p₀::State, λ₀::State=zero(q₀), μ₀::State=zero(q₀); kwargs...) = LDAE(ϑ, f, g, ḡ, ϕ, ψ, t₀, [q₀], [p₀], [λ₀], [μ₀]; kwargs...)
LDAE(ϑ, f, g, ḡ, ϕ, ψ, q₀::State, p₀::State, λ₀::State=zero(q₀), μ₀::State=zero(q₀); kwargs...) = LDAE(ϑ, f, g, ḡ, ϕ, ψ, 0.0, q₀, p₀, λ₀, μ₀; kwargs...)

const LDAEHT{HT,DT,TT,AT,ϑT,FT,GT,G̅T,ϕT,ψT,V̄T,F̄T,ΩT,∇T,PT} = LDAE{DT,TT,AT,ϑT,FT,GT,G̅T,ϕT,ψT,V̄T,F̄T,HT,ΩT,∇T,PT} # type alias for dispatch on Hamiltonian type parameter
const LDAE∇T{∇T,DT,TT,AT,ϑT,FT,GT,G̅T,ϕT,ψT,V̄T,F̄T,HT,ΩT,PT} = LDAE{DT,TT,AT,ϑT,FT,GT,G̅T,ϕT,ψT,V̄T,F̄T,HT,ΩT,∇T,PT} # type alias for dispatch on parameters type parameter
const LDAEΩT{ΩT,DT,TT,AT,ϑT,FT,GT,G̅T,ϕT,ψT,V̄T,F̄T,HT,∇T,PT} = LDAE{DT,TT,AT,ϑT,FT,GT,G̅T,ϕT,ψT,V̄T,F̄T,HT,ΩT,∇T,PT} # type alias for dispatch on symplectic two-form type parameter
const LDAEPT{PT,DT,TT,AT,ϑT,FT,GT,G̅T,ϕT,ψT,V̄T,F̄T,HT,ΩT,∇T} = LDAE{DT,TT,AT,ϑT,FT,GT,G̅T,ϕT,ψT,V̄T,F̄T,HT,ΩT,∇T,PT} # type alias for dispatch on parameters type parameter

Base.hash(dae::LDAE, h::UInt) = hash(dae.d, hash(dae.m,
          hash(dae.ϑ, hash(dae.f, hash(dae.g, hash(dae.ḡ, hash(dae.ϕ,
          hash(dae.ψ, hash(dae.v̄, hash(dae.f̄, hash(dae.h, hash(dae.Ω, hash(dae.∇H,
          hash(dae.t₀, hash(dae.q₀, hash(dae.p₀, hash(dae.λ₀, hash(dae.μ₀,
          hash(dae.parameters, hash(dae.periodicity, h))))))))))))))))))))

Base.:(==)(dae1::LDAE, dae2::LDAE) = (
                                dae1.d == dae2.d
                             && dae1.m == dae2.m
                             && dae1.ϑ == dae2.ϑ
                             && dae1.f == dae2.f
                             && dae1.g == dae2.g
                             && dae1.ḡ == dae2.ḡ
                             && dae1.ϕ == dae2.ϕ
                             && dae1.ψ == dae2.ψ
                             && dae1.v̄ == dae2.v̄
                             && dae1.f̄ == dae2.f̄
                             && dae1.h == dae2.h
                             && dae1.Ω == dae2.Ω
                             && dae1.∇H == dae2.∇H
                             && dae1.t₀ == dae2.t₀
                             && dae1.q₀ == dae2.q₀
                             && dae1.p₀ == dae2.p₀
                             && dae1.λ₀ == dae2.λ₀
                             && dae1.μ₀ == dae2.μ₀
                             && dae1.parameters == dae2.parameters
                             && dae1.periodicity == dae2.periodicity)

Base.similar(equ::LDAE, q₀, p₀, λ₀=get_λ₀(q₀, equ.λ₀), μ₀=get_λ₀(q₀, equ.μ₀); kwargs...) = similar(equ, equ.t₀, q₀, p₀, λ₀, μ₀; kwargs...)
Base.similar(equ::LDAE, t₀::Real, q₀::State, p₀::State, λ₀::State=get_λ₀(q₀, equ.λ₀), μ₀=get_λ₀(q₀, equ.μ₀); kwargs...) = similar(equ, t₀, [q₀], [p₀], [λ₀], [μ₀]; kwargs...)

function Base.similar(equ::LDAE, t₀::Real, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector;
                      v̄=equ.v̄, f̄=equ.f̄, h=equ.h, Ω=equ.Ω, ∇H=equ.∇H, parameters=equ.parameters, periodicity=equ.periodicity)
    @assert all([length(q) == ndims(equ) for q in q₀])
    @assert all([length(p) == ndims(equ) for p in p₀])
    @assert all([length(λ) == ndims(equ) for λ in λ₀])
    @assert all([length(μ) == ndims(equ) for μ in μ₀])
    LDAE(equ.ϑ, equ.f, equ.g, equ.ḡ, equ.ϕ, equ.ψ, t₀, q₀, p₀, λ₀, μ₀; v̄=v̄, f̄=f̄, h=h, Ω=Ω, ∇H=∇H, parameters=parameters, periodicity=periodicity)
end

@inline Base.ndims(equation::LDAE) = equation.d
@inline Base.axes(equation::LDAE) = axes(equation.q₀[begin])
@inline Common.nsamples(equ::LDAE) = length(eachindex(equ.q₀))
@inline Common.nconstraints(equation::LDAE) = equation.m
@inline Common.periodicity(equation::LDAE) = equation.periodicity

initial_conditions(equation::LDAE) = (equation.t₀, equation.q₀, equation.p₀, equation.λ₀, equation.μ₀)

hashamiltonian(::LDAEHT{<:Nothing}) = false
hashamiltonian(::LDAEHT{<:Function}) = true

hasgradientham(::LDAE∇T{<:Nothing}) = false
hasgradientham(::LDAE∇T{<:Function}) = true

hassymplecticform(::LDAEΩT{<:Nothing}) = false
hassymplecticform(::LDAEΩT{<:Function}) = true

hasparameters(::LDAEPT{<:Nothing}) = false
hasparameters(::LDAEPT{<:NamedTuple}) = true

_get_ϑ(equ::LDAE) = hasparameters(equ) ? (t,q,v,ϑ) -> equ.ϑ(t, q, v, ϑ, equ.parameters) : equ.ϑ
_get_f(equ::LDAE) = hasparameters(equ) ? (t,q,v,f) -> equ.f(t, q, v, f, equ.parameters) : equ.f
_get_g(equ::LDAE) = hasparameters(equ) ? (t,q,v,g) -> equ.g(t, q, v, g, equ.parameters) : equ.g
_get_ḡ(equ::LDAE) = hasparameters(equ) ? (t,q,λ,ḡ) -> equ.ḡ(t, q, λ, ḡ, equ.parameters) : equ.ḡ
_get_ϕ(equ::LDAE) = hasparameters(equ) ? (t,q,v,ϕ) -> equ.ϕ(t, q, v, ϕ, equ.parameters) : equ.ϕ
_get_ψ(equ::LDAE) = hasparameters(equ) ? (t,q,v,p,f,ψ) -> equ.ψ(t, q, v, p, f, ψ, equ.parameters) : equ.ψ
_get_v̄(equ::LDAE) = hasparameters(equ) ? (t,q,v) -> equ.v̄(t, q, v, equ.parameters) : equ.v̄
_get_f̄(equ::LDAE) = hasparameters(equ) ? (t,q,v,f) -> equ.f̄(t, q, v, f, equ.parameters) : equ.f̄
_get_h(equ::LDAE) = hasparameters(equ) ? (t,q) -> equ.h(t, q, equ.parameters) : equ.h
_get_∇(equ::LDAE) = hasparameters(equ) ? (t,q,∇H) -> equ.∇H(t, q, ∇H, equ.parameters) : equ.∇H
_get_Ω(equ::LDAE) = hasparameters(equ) ? (t,q,Ω) -> equ.Ω(t, q, Ω, equ.parameters) : equ.Ω


function get_function_tuple(equ::LDAE)
    names = (:ϑ, :f, :g, :ḡ, :ϕ, :ψ, :v̄, :f̄)
    equs  = (_get_ϑ(equ), _get_f(equ), _get_g(equ), _get_ḡ(equ), _get_ϕ(equ), _get_ψ(equ), _get_v̄(equ), _get_f̄(equ))

    if hashamiltonian(equ)
        names = (names..., :h)
        equs  = (equs..., _get_h(equ))
    end

    if hasgradientham(equ)
        names = (names..., :∇H)
        equs  = (equs..., _get_∇(equ))
    end

    if hassymplecticform(equ)
        names = (names..., :Ω)
        equs  = (equs..., _get_Ω(equ))
    end

    NamedTuple{names}(equs)
end
