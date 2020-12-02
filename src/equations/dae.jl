@doc raw"""
`DAE`: Differential Algebraic Equation

Defines a differential algebraic initial value problem
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t)) + u(t, q(t), \lambda(t)) , & q(t_{0}) &= q_{0} , \\
0 &= \phi (t, q(t), \lambda(t)) , & \lambda(t_{0}) &= \lambda_{0} ,
\end{aligned}
```
with vector field ``v``, projection ``u``, algebraic constraint ``\phi=0``,
initial conditions ``q_{0}`` and ``\lambda_{0}``, the dynamical variable ``q``
taking values in ``\mathbb{R}^{m}`` and the algebraic variable ``\lambda``
taking values in ``\mathbb{R}^{n}``.

### Fields

* `d`: dimension of dynamical variable ``q`` and the vector field ``v``
* `m`: dimension of algebraic variable ``\lambda`` and the constraint ``\phi``
* `n`: number of initial conditions
* `v`: function computing the vector field
* `u`: function computing the projection
* `ϕ`: algebraic constraint
* `v̄`: function computing an initial guess for the velocity field ``v``` (optional)
* `h`: function computing the Hamiltonian (optional)
* `t₀`: initial time
* `q₀`: initial condition for dynamical variable ``q``
* `λ₀`: initial condition for algebraic variable ``\lambda``

The function `v`, providing the vector field, takes three arguments,
`v(t, q, v)`, the functions `u` and `ϕ`, providing the projection and the
algebraic constraint take four arguments, `u(t, q, λ, u)` and `ϕ(t, q, λ, ϕ)`,
where `t` is the current time, `q` and `λ` are the current solution vectors,
and `v`, `u` and `ϕ` are the vectors which hold the result of evaluating the
vector field ``v``, the projection ``u`` and the algebraic constraint ``\phi``
on `t`, `q` and `λ`.

### Example

```julia
    function v(t, q, v)
        v[1] = q[1]
        v[2] = q[2]
    end

    function u(t, q, λ, u)
        u[1] = +λ[1]
        u[2] = -λ[1]
    end

    function ϕ(t, q, λ, ϕ)
        ϕ[1] = q[2] - q[1]
    end

    t₀ = 0.
    q₀ = [1., 1.]
    λ₀ = [0.]

    dae = DAE(v, u, ϕ, t₀, q₀, λ₀)

```
"""
struct DAE{dType <: Number, tType <: Number, vType <: Function, uType <: Function,
           ϕType <: Function, v̄Type <: Function, hType <: Union{Function,Nothing},
           pType <: Union{NamedTuple,Nothing}, N} <: AbstractEquationDAE{dType, tType}

    d::Int
    m::Int
    n::Int
    v::vType
    u::uType
    ϕ::ϕType
    v̄::v̄Type
    h::hType
    t₀::tType
    q₀::Array{dType, N}
    λ₀::Array{dType, N}
    parameters::pType
    periodicity::Vector{dType}

    function DAE(DT::DataType, N::Int, d::Int, m::Int, n::Int,
                 v::vType, u::uType, ϕ::ϕType, t₀::tType,
                 q₀::AbstractArray{dType}, λ₀::AbstractArray{dType};
                 v̄::v̄Type=v, h::hType=nothing, parameters::pType=nothing, periodicity=zeros(DT,d)) where {
                        dType <: Number, tType <: Number, vType <: Function,
                        uType <: Function, ϕType <: Function, v̄Type <: Function, 
                        hType <: Union{Function,Nothing}, pType <: Union{NamedTuple,Nothing}}

        @assert d == size(q₀,1)
        @assert m == size(λ₀,1)
        @assert n == size(q₀,2) == size(λ₀,2)
        @assert d ≥ m
        @assert ndims(q₀) == ndims(λ₀) == N ∈ (1,2)

        new{DT, tType, vType, uType, ϕType, v̄Type, hType, pType, N}(d, m, n, v, u, ϕ, v̄, h, t₀,
                convert(Array{DT}, q₀), convert(Array{DT}, λ₀), parameters, periodicity)
    end
end

function DAE(v, u, ϕ, t₀, q₀::AbstractArray{DT}, λ₀::AbstractArray{DT}; kwargs...) where {DT}
    DAE(DT, ndims(q₀), size(q₀,1), size(λ₀,1), size(q₀,2), v, u, ϕ, t₀, q₀, λ₀; kwargs...)
end

function DAE(v, u, ϕ, q₀, λ₀; kwargs...)
    DAE(v, u, ϕ, zero(eltype(q₀)), q₀, λ₀; kwargs...)
end

Base.hash(dae::DAE, h::UInt) = hash(dae.d, hash(dae.m, hash(dae.n, hash(dae.v,
        hash(dae.u, hash(dae.ϕ, hash(dae.v̄, hash(dae.h, hash(dae.t₀, hash(dae.q₀, hash(dae.λ₀,
        hash(dae.periodicity, hash(dae.parameters, h)))))))))))))

Base.:(==)(dae1::DAE, dae2::DAE) = (
                                dae1.d == dae2.d
                             && dae1.m == dae2.m
                             && dae1.n == dae2.n
                             && dae1.v == dae2.v
                             && dae1.u == dae2.u
                             && dae1.ϕ == dae2.ϕ
                             && dae1.v̄ == dae2.v̄
                             && dae1.h == dae2.h
                             && dae1.t₀ == dae2.t₀
                             && dae1.q₀ == dae2.q₀
                             && dae1.λ₀ == dae2.λ₀
                             && dae1.parameters == dae1.parameters
                             && dae1.periodicity == dae1.periodicity)

function Base.similar(dae::DAE, q₀, λ₀=get_λ₀(q₀, dae.λ₀); kwargs...)
    similar(dae, dae.t₀, q₀, λ₀; kwargs...)
end

function Base.similar(dae::DAE, t₀::TT, q₀::AbstractArray{DT}, λ₀::AbstractArray{DT}=get_λ₀(q₀, dae.λ₀);
                      v̄=dae.v̄, h=dae.h, parameters=dae.parameters, periodicity=dae.periodicity) where {DT  <: Number, TT <: Number}
    @assert dae.d == size(q₀,1)
    @assert dae.m == size(λ₀,1)
    DAE(dae.v, dae.u, dae.ϕ, t₀, q₀, λ₀; v̄=v̄, h=h, parameters=parameters, periodicity=periodicity)
end

@inline Base.ndims(equation::DAE) = equation.d
@inline CommonFunctions.nconstraints(equation::DAE) = equation.m
@inline CommonFunctions.periodicity(equation::DAE) = equation.periodicity

function get_function_tuple(equation::DAE{DT,TT,VT,UT,ϕT,V̄T,HT,Nothing}) where {DT, TT, VT, UT, ϕT, V̄T, HT}
    names = (:v,:u,:ϕ,:v̄)
    equs  = (equation.v, equation.u, equation.ϕ, equation.v̄)

    if HT != Nothing
        names = (names..., :h)
        equs  = (equs..., equation.h)
    end

    NamedTuple{names}(equs)
end

function get_function_tuple(equation::DAE{DT,TT,VT,UT,ϕT,V̄T,HT,PT}) where {DT, TT, VT, UT, ϕT, V̄T, HT, PT <: NamedTuple}
    vₚ = (t,q,v)   -> equation.v(t, q, v, equation.parameters)
    uₚ = (t,q,λ,u) -> equation.u(t, q, λ, u, equation.parameters)
    ϕₚ = (t,q,ϕ)   -> equation.ϕ(t, q, ϕ, equation.parameters)
    v̄ₚ = (t,q,v)   -> equation.v̄(t, q, v, equation.parameters)

    names = (:v, :u, :ϕ, :v̄)
    equs  = (vₚ, uₚ, ϕₚ, v̄ₚ)

    if HT != Nothing
        hₚ = (t,q,v) -> equation.h(t, q, v, equation.parameters)
        names = (names..., :h)
        equs  = (equs..., hₚ)
    end

    NamedTuple{names}(equs)
end
