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

### Parameters

* `DT <: Number`: data type
* `TT <: Real`: time step type
* `AT <: AbstractArray{DT}`: array type
* `vType <: Function`: type of `v`
* `fType <: Function`: type of `f`
* `ϕType <: Function`: type of `ϕ`
* `ψType <: Function`: type of `ψ`
* `hType <: OptionalFunction`: type of `h`
* `pType <: Union{NamedTuple,Nothing}`: parameters type

### Fields

* `d`: dimension of dynamical variables ``q`` and ``p`` as well as the vector fields ``v`` and ``f``
* `m`: dimension of algebraic variables ``\lambda`` and ``\gamma`` and the constraints ``\phi`` and ``\psi``
* `v`: tuple of functions computing the vector fields ``v_i``, ``i = 1 ... 3``
* `f`: tuple of functions computing the vector fields ``f_i``, ``i = 1 ... 3``
* `ϕ`: primary constraints
* `ψ`: secondary constraints
* `t₀`: initial time
* `q₀`: initial condition for dynamical variable ``q``
* `p₀`: initial condition for dynamical variable ``p``
* `λ₀`: initial condition for algebraic variable ``λ``
* `parameters`: either a `NamedTuple` containing the equations parameters or `nothing`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions

### Constructors

```julia
SPDAE(v, f, ϕ, ψ, t₀, q₀, p₀, λ₀; parameters=nothing, periodicity=zeros(DT,d))
SPDAE(v, f, ϕ, ψ, q₀::StateVector, p₀::StateVector, λ₀::StateVector; kwargs...) = SPDAE(v, f, ϕ, ψ, 0.0, q₀, p₀, λ₀; kwargs...)
SPDAE(v, f, ϕ, ψ, t₀, q₀::State, p₀::State, λ₀::State=zero(q₀); kwargs...) = SPDAE(v, f, ϕ, ψ, t₀, [q₀], [p₀], [λ₀]; kwargs...)
SPDAE(v, f, ϕ, ψ, q₀::State, p₀::State, λ₀::State=zero(q₀); kwargs...) = SPDAE(v, f, ϕ, ψ, 0.0, q₀, p₀, λ₀; kwargs...)
```

"""
struct SPDAE{dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
             vType <: Tuple, fType <: Tuple, ϕType <: Function, ψType <: Function,
             pType <: Union{NamedTuple,Nothing}} <: AbstractEquationPDAE{dType, tType}
    d::Int
    m::Int
    v::vType
    f::fType
    ϕ::ϕType
    ψ::ψType
    t₀::tType
    q₀::Vector{arrayType}
    p₀::Vector{arrayType}
    λ₀::Vector{arrayType}
    parameters::pType
    periodicity::arrayType

    function SPDAE(v::vType, f::fType, ϕ::ϕType, ψ::ψType,
                   t₀::tType, q₀::Vector{arrayType}, p₀::Vector{arrayType}, λ₀::Vector{arrayType};
                   parameters::pType=nothing, periodicity=zero(q₀[begin])) where {
                        dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
                        vType <: Tuple, fType <: Tuple, ϕType <: Function, ψType <: Function,
                        pType <: Union{NamedTuple,Nothing}}

        d = length(q₀[begin])
        m = length(λ₀[begin])

        @assert 2d ≥ m

        @assert length(q₀) == length(p₀) == length(λ₀)

        @assert all(length(q) == d for q in q₀)
        @assert all(length(p) == d for p in p₀)
        @assert all(length(λ) == m for λ in λ₀)

        @assert all([ndims(q) == ndims(p) == ndims(λ) for (q,p,λ) in zip(q₀,p₀,λ₀)])

        new{dType, tType, arrayType, vType, fType, ϕType, ψType, pType}(d, m, v, f, ϕ, ψ, t₀, q₀, p₀, λ₀, parameters, periodicity)
    end
end

SPDAE(v, f, ϕ, ψ, q₀::StateVector, p₀::StateVector, λ₀::StateVector; kwargs...) = SPDAE(v, f, ϕ, ψ, 0.0, q₀, p₀, λ₀; kwargs...)
SPDAE(v, f, ϕ, ψ, t₀, q₀::State, p₀::State, λ₀::State=zero(q₀); kwargs...) = SPDAE(v, f, ϕ, ψ, t₀, [q₀], [p₀], [λ₀]; kwargs...)
SPDAE(v, f, ϕ, ψ, q₀::State, p₀::State, λ₀::State=zero(q₀); kwargs...) = SPDAE(v, f, ϕ, ψ, 0.0, q₀, p₀, λ₀; kwargs...)

Base.hash(dae::SPDAE, h::UInt) = hash(dae.d, hash(dae.m,
                    hash(dae.v, hash(dae.f, hash(dae.ϕ, hash(dae.ψ,
                    hash(dae.t₀, hash(dae.q₀, hash(dae.p₀, h)))))))))

Base.:(==)(dae1::SPDAE, dae2::SPDAE) = (
                                dae1.d == dae2.d
                             && dae1.m == dae2.m
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
@inline Base.axes(equation::SPDAE) = axes(equation.q₀[begin])
@inline Common.nsamples(equation::SPDAE) = length(eachindex(equation.q₀))
@inline Common.nconstraints(equation::SPDAE) = equation.m
@inline Common.periodicity(equation::SPDAE) = equation.periodicity

function get_function_tuple(equation::SPDAE{DT,TT,AT,VT,FT,ϕT,ψT,Nothing}) where {DT, TT, AT, VT, FT, ϕT, ψT}
    NamedTuple{(:v, :f, :ϕ, :ψ)}((equation.v, equation.f, equation.ϕ, equation.ψ))
end

function get_function_tuple(equation::SPDAE{DT,TT,AT,VT,FT,ϕT,ψT,PT}) where {DT, TT, AT, VT, FT, ϕT, ψT, PT <: NamedTuple}
    vₚ = (t,q,p,v)   -> equation.v(t, q, p, v, equation.parameters)
    fₚ = (t,q,p,f)   -> equation.f(t, q, p, f, equation.parameters)
    ϕₚ = (t,q,p,ϕ)   -> equation.ϕ(t, q, p, ϕ, equation.parameters)
    ψₚ = (t,q,p,v,f,ψ) -> equation.ψ(t, q, p, v, f, ψ, equation.parameters)

    names = (:v, :f, :ϕ, :ψ)
    equs  = (vₚ, fₚ, ϕₚ, ψₚ)

    NamedTuple{names}(equs)
end
