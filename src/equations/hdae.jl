@doc raw"""
`HDAE`: Hamiltonian Differential Algebraic Equation *EXPERIMENTAL*

Defines a Hamiltonian differential algebraic initial value problem, that is
a canonical Hamiltonian system of equations subject to Dirac constraints,
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t), p(t)) + u(t, q(t), p(t), \lambda(t)) + \bar{g}(t, q(t), p(t), \lambda(t), \gamma(t)) , & q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), p(t)) + g(t, q(t), p(t), \lambda(t)) + \bar{f}(t, q(t), p(t), \lambda(t), \gamma(t)) , & p(t_{0}) &= p_{0} , \\
0 &= \phi (t, q(t), p(t)) , \\
0 &= \psi (t, q(t), p(t), \dot{q}(t), \dot{p}(t)) ,
\end{aligned}
```
with vector fields ``v``, ``u``, ``\bar{u}`` and ``f``, ``g``, ``\bar{g}``,
primary constraint ``\phi(q,p)=0`` and secondary constraint ``\psi(q,p,\lambda)=0``,
initial conditions ``(q_{0}, p_{0})``, the dynamical variables ``(q,p)``
taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}`` and
the algebraic variables ``(\lambda, \gamma)`` taking values in
``\mathbb{R}^{n} \times \mathbb{R}^{d}``.

### Parameters

* `DT <: Number`: data type
* `TT <: Real`: time step type
* `AT <: AbstractArray{DT}`: array type
* `vType <: Function`: type of `v`
* `fType <: Function`: type of `f`
* `uType <: Function`: type of `u`
* `gType <: Function`: type of `g`
* `ϕType <: Function`: type of `ϕ`
* `ūType <: Function`: type of `ū`
* `ḡType <: Function`: type of `ḡ`
* `ψType <: Function`: type of `ψ`
* `PType <: Function`: type of `P`
* `v̄Type <: Function`: type of `v̄`
* `f̄Type <: Function`: type of `f̄`
* `hamType <: Function`: Hamiltonian type
* `invType <: OptionalNamedTuple`: invariants type
* `parType <: OptionalNamedTuple`: parameters type
* `perType <: OptionalArray{AT}`: periodicity type

### Fields

* `d`: dimension of dynamical variables ``q`` and ``p`` as well as the vector fields ``v`` and ``f``
* `m`: dimension of algebraic variables ``\lambda`` and ``\gamma`` and the constraints ``\phi`` and ``\psi``
* `v`: function computing the Hamiltonian vector field ``v``
* `f`: function computing the Hamiltonian vector field ``f``
* `u`: function computing the projection for ``q``
* `g`: function computing the primary projection field ``g``
* `ϕ`: primary constraints
* `ū`: function computing the secondary projection field ``\bar{u}`` (optional)
* `ḡ`: function computing the secondary projection field ``\bar{g}`` (optional)
* `ψ`: secondary constraints (optional)
* `P`: function computing the Poisson matrix ``P``
* `v̄`: function computing an initial guess for the velocity field ``v``` (optional, defaults to `v`)
* `f̄`: function computing an initial guess for the force field ``f`` (optional, defaults to `f`)
* `t₀`: initial time (optional)
* `q₀`: initial condition for dynamical variable ``q``
* `p₀`: initial condition for dynamical variable ``p``
* `λ₀`: initial condition for algebraic variable ``λ``
* `hamiltonian`: function computing the Hamiltonian ``H``
* `invariants`: either a `NamedTuple` containing the equation's invariants or `nothing`
* `parameters`: either a `NamedTuple` containing the equation's parameters or `nothing`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions

### Constructors

```julia
HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, P, t₀, q₀, p₀, λ₀, hamiltonian, invariants, parameters, periodicity)

HDAE(v, f, u, g, ϕ, h, t₀, q₀::StateVector, p₀::StateVector, λ₀::StateVector; kwargs...)
HDAE(v, f, u, g, ϕ, h, q₀::StateVector, p₀::StateVector, λ₀::StateVector; kwargs...)
HDAE(v, f, u, g, ϕ, h, t₀, q₀::State, p₀::State, λ₀::State; kwargs...)
HDAE(v, f, u, g, ϕ, h, q₀::State, p₀::State, λ₀::State; kwargs...)

HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, h, t₀, q₀::StateVector, p₀::StateVector, λ₀::StateVector; kwargs...)
HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, h, q₀::StateVector, p₀::StateVector, λ₀::StateVector; kwargs...)
HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, h, t₀, q₀::State, p₀::State, λ₀::State; kwargs...)
HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, h, q₀::State, p₀::State, λ₀::State; kwargs...)
```

"""
struct HDAE{dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
            vType <: Function, fType <: Function,
            uType <: Function, gType <: Function, ϕType <: Function,
            ūType <: OptionalFunction, ḡType <: OptionalFunction, ψType <: OptionalFunction,
            v̄Type <: Function, f̄Type <: Function,
            PType <: Function,
            hamType <: Function,
            invType <: OptionalNamedTuple,
            parType <: OptionalNamedTuple,
            perType <: OptionalArray{arrayType}} <: AbstractEquationPDAE{dType, tType}

    d::Int
    m::Int

    v::vType
    f::fType
    u::uType
    g::gType
    ϕ::ϕType
    ū::ūType
    ḡ::ḡType
    ψ::ψType
    v̄::v̄Type
    f̄::f̄Type
    P::PType

    t₀::tType
    q₀::Vector{arrayType}
    p₀::Vector{arrayType}
    λ₀::Vector{arrayType}
    μ₀::Vector{arrayType}

    hamiltonian::hamType
    invariants::invType
    parameters::parType
    periodicity::perType

    function HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, P,
                  t₀::tType, q₀::Vector{arrayType}, p₀::Vector{arrayType}, λ₀::Vector{arrayType}, μ₀::Vector{arrayType},
                  hamiltonian, invariants, parameters, periodicity) where {
                        dType <: Number, tType <: Real, arrayType <: AbstractArray{dType}}

        d = length(q₀[begin])
        m = length(λ₀[begin])

        @assert 2d ≥ m

        @assert length(q₀) == length(p₀) == length(λ₀)

        @assert all(length(q) == d for q in q₀)
        @assert all(length(p) == d for p in p₀)
        @assert all(length(λ) == m for λ in λ₀)
        @assert all(length(μ) == m for μ in μ₀)

        @assert all([ndims(q) == ndims(p) == ndims(λ) == ndims(μ) for (q,p,λ,μ) in zip(q₀,p₀,λ₀,μ₀)])

        new{dType, tType, arrayType,
            typeof(v), typeof(f),
            typeof(u), typeof(g), typeof(ϕ),
            typeof(ū), typeof(ḡ), typeof(ψ),
            typeof(v̄), typeof(f̄), typeof(P),
            typeof(hamiltonian), typeof(invariants), typeof(parameters), typeof(periodicity)}(
                d, m, v, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, P, t₀, q₀, p₀, λ₀, μ₀,
                hamiltonian, invariants, parameters, periodicity)
    end
end

_HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, hamiltonian, t₀, q₀, p₀, λ₀, μ₀; v̄=v, f̄=f, P=symplectic_matrix, invariants=nothing, parameters=nothing, periodicity=nothing) = HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, P, t₀, q₀, p₀, λ₀, μ₀, hamiltonian, invariants, parameters, periodicity)

HDAE(v, f, u, g, ϕ, h, t₀, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...) = _HDAE(v, f, u, g, ϕ, nothing, nothing, nothing, h, t₀, q₀, p₀, λ₀, μ₀; kwargs...)
HDAE(v, f, u, g, ϕ, h, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...) = HDAE(v, f, u, g, ϕ, h, 0.0, q₀, p₀, λ₀, μ₀; kwargs...)
HDAE(v, f, u, g, ϕ, h, t₀, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...) = HDAE(v, f, u, g, ϕ, h, t₀, [q₀], [p₀], [λ₀], [μ₀]; kwargs...)
HDAE(v, f, u, g, ϕ, h, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...) = HDAE(v, f, u, g, ϕ, h, 0.0, q₀, p₀, λ₀, μ₀; kwargs...)

HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, h, t₀, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...) = _HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, h, t₀, q₀, p₀, λ₀, μ₀; kwargs...)
HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, h, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...) = HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, h, 0.0, q₀, p₀, λ₀, μ₀; kwargs...)
HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, h, t₀, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...) = HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, h, t₀, [q₀], [p₀], [λ₀], [μ₀]; kwargs...)
HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, h, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...) = HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, h, 0.0, q₀, p₀, λ₀, μ₀; kwargs...)

const HDAEpsiType{psiT,DT,TT,AT,VT,FT,UT,GT,ΦT,ŪT,ḠT,V̄T,F̄T,PT,hamT,invT,parT,perT} = HDAE{DT,TT,AT,VT,FT,UT,GT,ΦT,ŪT,ḠT,psiT,V̄T,F̄T,PT,hamT,invT,parT,perT} # type alias for dispatch on secondary constraint type parameter
const HDAEinvType{invT,DT,TT,AT,VT,FT,UT,GT,ΦT,ŪT,ḠT,psiT,V̄T,F̄T,PT,hamT,parT,perT} = HDAE{DT,TT,AT,VT,FT,UT,GT,ΦT,ŪT,ḠT,psiT,V̄T,F̄T,PT,hamT,invT,parT,perT} # type alias for dispatch on invariants type parameter
const HDAEparType{parT,DT,TT,AT,VT,FT,UT,GT,ΦT,ŪT,ḠT,psiT,V̄T,F̄T,PT,hamT,invT,perT} = HDAE{DT,TT,AT,VT,FT,UT,GT,ΦT,ŪT,ḠT,psiT,V̄T,F̄T,PT,hamT,invT,parT,perT} # type alias for dispatch on parameters type parameter
const HDAEperType{perT,DT,TT,AT,VT,FT,UT,GT,ΦT,ŪT,ḠT,psiT,V̄T,F̄T,PT,hamT,invT,parT} = HDAE{DT,TT,AT,VT,FT,UT,GT,ΦT,ŪT,ḠT,psiT,V̄T,F̄T,PT,hamT,invT,parT,perT} # type alias for dispatch on periodicity type parameter

Base.hash(dae::HDAE, h::UInt) = hash(dae.d, hash(dae.m,
                        hash(dae.v, hash(dae.f,
                        hash(dae.u, hash(dae.g, hash(dae.ϕ,
                        hash(dae.ū, hash(dae.ḡ, hash(dae.ψ,
                        hash(dae.v̄, hash(dae.f̄, hash(dae.P,
                        hash(dae.t₀, hash(dae.q₀, hash(dae.p₀, hash(dae.λ₀, hash(dae.μ₀,
                        hash(dae.hamiltonian, hash(dae.invariants, hash(dae.parameters, hash(dae.periodicity, h))))))))))))))))))))))

Base.:(==)(dae1::HDAE, dae2::HDAE) = (
                                dae1.d == dae2.d
                             && dae1.m == dae2.m
                             && dae1.v == dae2.v
                             && dae1.f == dae2.f
                             && dae1.u == dae2.u
                             && dae1.g == dae2.g
                             && dae1.ϕ == dae2.ϕ
                             && dae1.ū == dae2.ū
                             && dae1.ḡ == dae2.ḡ
                             && dae1.ψ == dae2.ψ
                             && dae1.v̄ == dae2.v̄
                             && dae1.f̄ == dae2.f̄
                             && dae1.P == dae2.P
                             && dae1.t₀ == dae2.t₀
                             && dae1.q₀ == dae2.q₀
                             && dae1.p₀ == dae2.p₀
                             && dae1.λ₀ == dae2.λ₀
                             && dae1.μ₀ == dae2.μ₀
                             && dae1.hamiltonian == dae2.hamiltonian
                             && dae1.invariants  == dae2.invariants
                             && dae1.parameters  == dae2.parameters
                             && dae1.periodicity == dae2.periodicity)

function Base.similar(equ::HDAE, t₀::Real, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector; parameters=equ.parameters)
    @assert all([length(q) == equ.d for q in q₀])
    @assert all([length(p) == equ.d for p in p₀])
    @assert all([length(λ) == equ.m for λ in λ₀])
    @assert all([length(μ) == equ.m for μ in μ₀])
    _HDAE(equ.v, equ.f, equ.u, equ.g, equ.ϕ, equ.ū, equ.ḡ, equ.ψ, equ.hamiltonian, t₀, q₀, p₀, λ₀, μ₀;
         v̄=equ.v̄, f̄=equ.f̄, P=equ.P, invariants=equ.invariants, parameters=parameters, periodicity=equ.periodicity)
end

Base.similar(equ::HDAE, q₀, p₀, λ₀=get_λ₀(q₀, equ.λ₀), μ₀=get_λ₀(λ₀, equ.μ₀); kwargs...) = similar(equ, equ.t₀, q₀, p₀, λ₀, μ₀; kwargs...)
Base.similar(equ::HDAE, t₀::Real, q₀::State, p₀::State, λ₀::State=get_λ₀(q₀, equ.λ₀), μ₀::State=get_λ₀(λ₀, equ.μ₀); kwargs...) = similar(equ, t₀, [q₀], [p₀], [λ₀], [μ₀]; kwargs...)

hashamiltonian(::HDAE) = true

hassecondary(::HDAEpsiType{<:Nothing}) = false
hassecondary(::HDAEpsiType{<:Function}) = true

hasinvariants(::HDAEinvType{<:Nothing}) = false
hasinvariants(::HDAEinvType{<:NamedTuple}) = true

hasparameters(::HDAEparType{<:Nothing}) = false
hasparameters(::HDAEparType{<:NamedTuple}) = true

hasperiodicity(::HDAEperType{<:Nothing}) = false
hasperiodicity(::HDAEperType{<:AbstractArray}) = true

@inline Base.axes(equation::HDAE) = axes(equation.q₀[begin])
@inline Base.ndims(equation::HDAE) = equation.d
@inline Common.nsamples(equation::HDAE) = length(eachindex(equation.q₀))
@inline Common.nconstraints(equation::HDAE) = equation.m

@inline Common.periodicity(equation::HDAE) = hasperiodicity(equation) ? equation.periodicity : zero(equation.q₀[begin])
@inline initial_conditions(equation::HDAE) = (equation.t₀, equation.q₀, equation.p₀, equation.λ₀)

_get_v(equ::HDAE)  = hasparameters(equ) ? (t,q,p,v)     -> equ.v(t, q, p, v, equ.parameters) : equ.v
_get_f(equ::HDAE)  = hasparameters(equ) ? (t,q,p,f)     -> equ.f(t, q, p, f, equ.parameters) : equ.f
_get_u(equ::HDAE)  = hasparameters(equ) ? (t,q,p,λ,u)   -> equ.u(t, q, p, λ, u, equ.parameters) : equ.u
_get_g(equ::HDAE)  = hasparameters(equ) ? (t,q,p,λ,g)   -> equ.g(t, q, p, λ, g, equ.parameters) : equ.g
_get_ϕ(equ::HDAE)  = hasparameters(equ) ? (t,q,p,ϕ)     -> equ.ϕ(t, q, p, ϕ, equ.parameters) : equ.ϕ
_get_ū(equ::HDAE)  = hasparameters(equ) ? (t,q,p,λ,ū)   -> equ.ū(t, q, p, λ, ū, equ.parameters) : equ.ū
_get_ḡ(equ::HDAE)  = hasparameters(equ) ? (t,q,p,λ,ḡ)   -> equ.ḡ(t, q, p, λ, ḡ, equ.parameters) : equ.ḡ
_get_ψ(equ::HDAE)  = hasparameters(equ) ? (t,q,p,v,f,ψ) -> equ.ψ(t, q, p, v, f, ψ, equ.parameters) : equ.ψ
_get_v̄(equ::HDAE)  = hasparameters(equ) ? (t,q,p,v)     -> equ.v̄(t, q, p, v, equ.parameters) : equ.v̄
_get_f̄(equ::HDAE)  = hasparameters(equ) ? (t,q,p,f)     -> equ.f̄(t, q, p, f, equ.parameters) : equ.f̄
_get_h(equ::HDAE)  = hasparameters(equ) ? (t,q,p)       -> equ.hamiltonian(t, q, p, equ.parameters) : equ.hamiltonian
_get_P(equ::HDAE)  = hasparameters(equ) ? (t,q,p,P)     -> equ.P(t, q, p, P, equ.parameters) : equ.P

function get_function_tuple(equ::HDAE)
    names = (:v, :f, :u, :g, :ϕ,)
    equs  = (_get_v(equ), _get_f(equ), 
             _get_u(equ), _get_g(equ), _get_ϕ(equ))

    if hassecondary(equ)
        names = (names..., :ū, :ḡ, :ψ)
        equs  = (equs..., _get_ū(equ), _get_ḡ(equ), _get_ψ(equ))
    end

    names = (names..., :h, :v̄, :f̄)
    equs  = (equs..., _get_h(equ), _get_v̄(equ), _get_f̄(equ))

    NamedTuple{names}(equs)
end
