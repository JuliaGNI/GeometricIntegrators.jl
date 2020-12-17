@doc raw"""
`SPSDE`: Stratonovich Split Partitioned Stochastic Differential Equation

Defines a partitioned stochastic differential initial value problem
```math
\begin{aligned}
dq (t) &=   v(t, q(t)) \, dt + B(t, q(t)) \circ dW , & q(t_{0}) &= q_{0} , \\
dp (t) &= [ f_1(t, q(t)) + f_2(t, q(t)) ] \, dt + [ G_1(t, q(t)) + G_2(t, q(t)) ] \circ dW , & p(t_{0}) &= p_{0} ,
\end{aligned}
```
with the drift vector fields ``v`` and ``f_i``, diffusion matrices ``B`` and ``G_i``,
initial conditions ``q_{0}`` and ``p_{0}``, the dynamical variables ``(q,p)`` taking
values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}``, and the m-dimensional Wiener process W

### Fields

* `d`:  dimension of dynamical variable ``q`` and the vector fields ``vi``
* `m`:  dimension of the Wiener process
* `ni`: number of initial conditions
* `ns`: number of sample paths
* `v` :  function computing the drift vector field for the position variable ``q``
* `f1`:  function computing the drift vector field for the momentum variable ``p``
* `f2`:  function computing the drift vector field for the momentum variable ``p``
* `B` :  function computing the d x m diffusion matrix for the position variable ``q``
* `G1`:  function computing the d x m diffusion matrix for the momentum variable ``p``
* `G2`:  function computing the d x m diffusion matrix for the momentum variable ``p``
* `t₀`: initial time
* `q₀`: initial condition for dynamical variable ``q`` (may be a random variable itself)
* `p₀`: initial condition for dynamical variable ``p`` (may be a random variable itself)


The functions `v`, `f`, `B` and `G`, providing the drift vector fields and diffusion matrices, take four arguments,
`v(t, q, p, v)`, `f(t, q, p, f)`, `B(t, q, p,  B)` and `G(t, q, p, G)`, where `t` is the current time, `(q, p)` is the
current solution vector, and `v`, `f`, `B` and `G` are the variables which hold the result
of evaluating the vector fields ``v``, ``f`` and the matrices ``B``, ``G`` on `t` and `(q,p)`.

### Example

```julia
    function v(λ, t, q, v)
        v[1] = λ*q[1]
        v[2] = λ*q[2]
    end

    function B(μ, t, q, B)
        B[1] = μ*q[1]
        B[2] = μ*q[2]
    end

    t₀ = 0.
    q₀ = [1., 1.]
    λ  = 2.
    μ  = 1.

    v_sde = (t, q, v) -> v(λ, t, q, v)
    B_sde = (t, q, B) -> B(μ, t, q, B)

    sde = SDE(v_sde, B_sde, t₀, q₀)
```
"""
struct SPSDE{dType <: Number, tType <: Real, vType <: Function, f1Type <: Function, f2Type <: Function, BType <: Function, G1Type <: Function, G2Type <: Function, pType, N} <: AbstractEquationPSDE{dType, tType}
    d::Int
    m::Int
    ni::Int
    ns::Int
    v ::vType
    f1::f1Type
    f2::f2Type
    B ::BType
    G1::G1Type
    G2::G2Type
    t₀::tType
    q₀::Array{dType, N}           #Initial condition: N=1 - single deterministic, N=2 - single random or multiple deterministic, N=3 - multiple deterministic
    p₀::Array{dType, N}
    parameters::pType
    periodicity::Vector{dType}

    function SPSDE(m, ns, v::vType, f1::f1Type, f2::f2Type,
                    B::BType, G1::G1Type, G2::G2Type, t₀::tType,
                    q₀::AbstractArray{dType,N}, p₀::AbstractArray{dType,N};
                    parameters=nothing, periodicity=zeros(dType,size(q₀,1))) where {
                    dType <: Number, tType <: Real, vType <: Function, f1Type <: Function, f2Type <: Function,
                    BType <: Function, G1Type <: Function, G2Type <: Function, N}

        @assert size(q₀)  == size(p₀)

        d  = size(q₀,1)
        ni = size(q₀,2)

        @assert N ∈ (1,2)
        @assert ni ≥ 1
        @assert ns ≥ 1
        @assert ni == 1 || ns == 1
        # either multiple deterministic initial conditions and one sample path
        # or one deterministic initial condition and multiple sample paths

        if !(length(periodicity) == d)
            periodicity = zeros(dType, d)
        end

        new{dType,tType,vType,f1Type,f2Type,BType,G1Type,G2Type,typeof(parameters),N}(d, m, ni, ns, v, f1, f2, B, G1, G2, t₀, q₀, p₀, parameters, periodicity)
    end
end


function SPSDE(m::Int, ns::Int, v::Function, f1::Function, f2::Function,
                                B::Function, G1::Function, G2::Function,
                                q₀::AbstractArray{DT}, p₀::AbstractArray{DT}; kwargs...) where {DT}
    SPSDE(m, ns, v, f1, f2, B, G1, G2, zero(DT), q₀, p₀; kwargs...)
end


Base.hash(sde::SPSDE, h::UInt) = hash(sde.d, hash(sde.m, hash(sde.ni, hash(sde.ns,
                                 hash(sde.v, hash(sde.f1, hash(sde.f2, hash(sde.B, hash(sde.G1, hash(sde.G2,
                                 hash(sde.t₀, hash(sde.q₀, hash(sde.p₀, hash(sde.parameters, hash(sde.periodicity, h)))))))))))))))

Base.:(==)(sde1::SPSDE, sde2::SPSDE) = (
                                sde1.d == sde2.d
                             && sde1.m == sde2.m
                             && sde1.ni == sde2.ni
                             && sde1.ns == sde2.ns
                             && sde1.v  == sde2.v
                             && sde1.f1 == sde2.f1
                             && sde1.f2 == sde2.f2
                             && sde1.B  == sde2.B
                             && sde1.G1 == sde2.G1
                             && sde1.G2 == sde2.G2
                             && sde1.t₀ == sde2.t₀
                             && sde1.q₀ == sde2.q₀
                             && sde1.p₀ == sde2.p₀
                             && sde1.parameters  == sde2.parameters
                             && sde1.periodicity == sde2.periodicity)

function Base.similar(sde::SPSDE, q₀::AbstractArray, p₀::AbstractArray, ns::Int)
    similar(sde, sde.t₀, q₀, p₀, ns)
end

function Base.similar(sde::SPSDE, q₀::AbstractArray, p₀::AbstractArray)
    similar(sde, sde.t₀, q₀, p₀)
end

function Base.similar(sde::SPSDE, t₀::TT, q₀::AbstractArray{DT,N}, p₀::AbstractArray{DT,N}, ns::Int=(N > 1 ? 1 : sde.ns)) where {DT <: Number, TT <: Number, N}
    @assert size(q₀) == size(p₀)
    @assert sde.d == size(q₀,1)
    SPSDE(sde.m, ns, sde.v, sde.f1, sde.f2, sde.B, sde.G1, sde.G2, t₀, q₀, p₀; parameters=sde.parameters, periodicity=sde.periodicity)
end

@inline Base.ndims(sde::SPSDE) = sde.d

@inline Common.periodicity(equation::SPSDE) = equation.periodicity

function get_function_tuple(equation::SPSDE{DT,TT,VT,F1T,F2T,BT,G1T,G2T,Nothing}) where {DT, TT, VT, F1T, F2T, BT, G1T, G2T}
    NamedTuple{(:v,:f1,:f2,:B,:G1,:G2)}((equation.v, equation.f1, equation.f2, equation.B, equation.G1, equation.G2))
end

function get_function_tuple(equation::SPSDE{DT,TT,VT,F1T,F2T,BT,G1T,G2T,PT}) where {DT, TT, VT, F1T, F2T, BT, G1T, G2T, PT <: NamedTuple}
    vₚ  = (t,q,p,v) -> equation.v(t, q, p, v, equation.parameters)
    f1ₚ = (t,q,p,f) -> equation.f1(t, q, p, f, equation.parameters)
    f2ₚ = (t,q,p,f) -> equation.f2(t, q, p, f, equation.parameters)
    Bₚ  = (t,q,p,B) -> equation.B(t, q, p, B, equation.parameters)
    G1ₚ = (t,q,p,G) -> equation.G1(t, q, p, G, equation.parameters)
    G2ₚ = (t,q,p,G) -> equation.G2(t, q, p, G, equation.parameters)

    names = (:v,:f1,:f2,:B,:G1,:G2)
    equs  = (vₚ, f1ₚ, f2ₚ, Bₚ, G1ₚ, G2ₚ)

    NamedTuple{names}(equs)
end
