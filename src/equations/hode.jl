@doc raw"""
`HODE`: Hamiltonian Ordinary Differential Equation *EXPERIMENTAL*

Defines a Hamiltonian ordinary differential initial value problem, that is
a canonical Hamiltonian system of equations,
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t), p(t)) , & q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), p(t)) , & p(t_{0}) &= p_{0} ,
\end{aligned}
```
with vector fields ``v`` and ``f``, given by
```math
\begin{aligned}
v &=   \frac{\partial H}{\partial p} , &
f &= - \frac{\partial H}{\partial q} ,
\end{aligned}
```
initial conditions ``(q_{0}, p_{0})`` and the dynamical variables ``(q,p)``
taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}``.

### Fields

* `d`: dimension of dynamical variables ``q`` and ``p`` as well as the vector fields ``v`` and ``f``
* `n`: number of initial conditions
* `v`: function computing the vector field ``v``
* `f`: function computing the vector field ``f``
* `h`: function computing the Hamiltonian ``H``
* `t₀`: initial time
* `q₀`: initial condition for dynamical variable ``q``
* `p₀`: initial condition for dynamical variable ``p``

"""
struct HODE{dType <: Number, tType <: Number,
            vType <: Function, fType <: Function, hType <: Function,
            pType <: Union{NamedTuple,Nothing}, N} <: AbstractEquationPODE{dType, tType}

    d::Int
    n::Int
    v::vType
    f::fType
    h::hType
    t₀::tType
    q₀::Array{dType, N}
    p₀::Array{dType, N}
    parameters::pType
    periodicity::Vector{dType}

    function HODE(DT::DataType, N::Int, d::Int, n::Int, v::vType, f::fType, h::hType,
                  t₀::tType, q₀::AbstractArray{dType}, p₀::AbstractArray{dType};
                  parameters::pType=nothing, periodicity=zeros(DT,d)) where {
                        dType <: Number, tType <: Number, vType <: Function,
                        fType <: Function, hType <: Function,
                        pType <: Union{NamedTuple,Nothing}}

        @assert d == size(q₀,1) == size(p₀,1)
        @assert n == size(q₀,2) == size(p₀,2)
        @assert ndims(q₀) == ndims(p₀) == N ∈ (1,2)

        new{DT, tType, vType, fType, hType, pType, N}(d, n, v, f, h, t₀,
                convert(Array{DT}, q₀), convert(Array{DT}, p₀),
                parameters, periodicity)
    end
end

function HODE(v, f, h, t₀, q₀::AbstractArray{DT}, p₀::AbstractArray{DT}; kwargs...) where {DT}
    HODE(DT, ndims(q₀), size(q₀,1), size(q₀,2), v, f, h, t₀, q₀, p₀; kwargs...)
end

function HODE(v, f, h, q₀, p₀; kwargs...)
    HODE(v, f, h, zero(eltype(q₀)), q₀, p₀; kwargs...)
end

Base.hash(ode::HODE, h::UInt) = hash(ode.d, hash(ode.n, hash(ode.v, hash(ode.f, hash(ode.h,
        hash(ode.t₀, hash(ode.q₀, hash(ode.p₀, hash(ode.periodicity, hash(ode.parameters, h))))))))))

Base.:(==)(ode1::HODE, ode2::HODE) = (
                                ode1.d == ode2.d
                             && ode1.n == ode2.n
                             && ode1.v == ode2.v
                             && ode1.f == ode2.f
                             && ode1.h == ode2.h
                             && ode1.t₀ == ode2.t₀
                             && ode1.q₀ == ode2.q₀
                             && ode1.p₀ == ode2.p₀
                             && ode1.parameters == ode2.parameters
                             && ode1.periodicity == ode2.periodicity)

function Base.similar(ode::HODE, q₀, p₀; kwargs...)
    similar(ode, ode.t₀, q₀, p₀; kwargs...)
end

function Base.similar(ode::HODE, t₀::TT, q₀::AbstractArray{DT}, p₀::AbstractArray{DT};
                      parameters=ode.parameters, periodicity=ode.periodicity) where {DT  <: Number, TT <: Number}
    @assert ode.d == size(q₀,1) == size(p₀,1)
    HODE(ode.v, ode.f, ode.h, t₀, q₀, p₀; parameters=parameters, periodicity=periodicity)
end

@inline Base.ndims(ode::HODE) = ode.d

@inline CommonFunctions.periodicity(equation::HODE) = equation.periodicity

function get_function_tuple(equation::HODE{DT,TT,VT,FT,HT,Nothing}) where {DT, TT, VT, FT, HT}
    NamedTuple{(:v,:f,:h)}((equation.v, equation.f, equation.h))
end

function get_function_tuple(equation::HODE{DT,TT,VT,FT,HT,PT}) where {DT, TT, VT, FT, HT, PT <: NamedTuple}
    vₚ = (t,q,p,v) -> equation.v(t, q, p, v, equation.parameters)
    fₚ = (t,q,p,f) -> equation.f(t, q, p, f, equation.parameters)
    hₚ = (t,q,p) -> equation.h(t, q, p, equation.parameters)

    names = (:v, :f, :h)
    equs  = (vₚ, fₚ, hₚ)

    NamedTuple{names}(equs)
end
