@doc raw"""
`PSDE`: Stratonovich Partitioned Stochastic Differential Equation

Defines a partitioned stochastic differential initial value problem
```math
\begin{align*}
\dq (t) &= v(t, q(t)) \, dt + B(t, q(t)) \circ dW , & q(t_{0}) &= q_{0} ,
\dp (t) &= f(t, q(t)) \, dt + G(t, q(t)) \circ dW , & p(t_{0}) &= p_{0}
\end{align*}
```
with the drift vector fields ``v`` and ``f``, diffusion matrices ``B`` and ``G``,
initial conditions ``q_{0}`` and ``p_{0}``, the dynamical variables ``(q,p)`` taking
values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}``, and the m-dimensional Wiener process W

### Fields

* `d`:  dimension of dynamical variable ``q`` and the vector field ``v``
* `m`:  dimension of the Wiener process
* `ni`: number of initial conditions
* `ns`: number of sample paths
* `v`:  function computing the drift vector field for the position variable q
* `f`:  function computing the drift vector field for the momentum variable p
* `B`:  function computing the d x m diffusion matrix for the position variable q
* `G`:  function computing the d x m diffusion matrix for the momentum variable p
* `t₀`: initial time
* `q₀`: initial condition for dynamical variable ``q`` (may be a random variable itself)
* `p₀`: initial condition for dynamical variable ``p`` (may be a random variable itself)


The functions `v`, `f`, 'B' and `G`, providing the drift vector fields and diffusion matrices, take four arguments,
`v(t, q, p, v)`, `f(t, q, p, f)`, `B(t, q, p,  B)` and `G(t, q, p, G)`, where `t` is the current time, `(q, p)` is the
current solution vector, and `v`, `f`, 'B' and `G` are the variables which hold the result
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
struct PSDE{dType <: Number, tType <: Number, vType <: Function, fType <: Function,
            BType <: Function, GType <: Function, pType <: Union{Tuple,Nothing}, N} <: AbstractEquationPSDE{dType, tType}

    d::Int
    m::Int
    ni::Int
    ns::Int
    v::vType
    f::fType
    B::BType
    G::GType
    t₀::tType
    q₀::Array{dType, N}           #Initial condition: N=1 - single deterministic, N=2 - single random or multiple deterministic, N=3 - multiple deterministic
    p₀::Array{dType, N}
    parameters::pType
    periodicity::Vector{dType}

    function PSDE(m, ns, v::vType, f::fType, B::BType, G::GType,
                  t₀::tType, q₀::DenseArray{dType}, p₀::DenseArray{dType};
                  parameters=nothing, periodicity=zeros(dType,size(q₀,1))) where {
                        dType <: Number, tType <: Number,
                        vType <: Function, fType <: Function,
                        BType <: Function, GType <: Function}

        @assert size(q₀)  == size(p₀)
        @assert ndims(q₀) == ndims(p₀)

        N  = ndims(q₀)
        d  = size(q₀,1)
        ni = size(q₀,2)

        @assert N ∈ (1,2)
        @assert ni ≥ 1
        @assert ns ≥ 1
        @assert ni == 1 || ns == 1
        # either multiple deterministic initial conditions and one sample path
        # or one deterministic initial condition and multiple sample paths

        new{dType,tType,vType,fType,BType,GType,typeof(parameters),N}(d, m, ni, ns, v, f, B, G, t₀, q₀, p₀, parameters, periodicity)
    end
end


function PSDE(m::Int, ns::Int, v::Function, f::Function, B::Function, G::Function, q₀::DenseArray{DT}, p₀::DenseArray{DT}; kwargs...) where {DT}
    PSDE(m, ns, v, f, B, G, zero(DT), q₀, p₀; kwargs...)
end


Base.hash(sde::PSDE, h::UInt) = hash(sde.d, hash(sde.m, hash(sde.ni, hash(sde.ns,
                                hash(sde.v, hash(sde.f, hash(sde.B, hash(sde.G,
                                hash(sde.t₀, hash(sde.q₀, hash(sde.p₀,
                                hash(sde.parameters, hash(sde.periodicity, h)))))))))))))

Base.:(==)(sde1::PSDE, sde2::PSDE) = (
                                sde1.d == sde2.d
                             && sde1.m == sde2.m
                             && sde1.ni == sde2.ni
                             && sde1.ns == sde2.ns
                             && sde1.v == sde2.v
                             && sde1.f == sde2.f
                             && sde1.B == sde2.B
                             && sde1.G == sde2.G
                             && sde1.t₀ == sde2.t₀
                             && sde1.q₀ == sde2.q₀
                             && sde1.p₀ == sde2.p₀
                             && sde1.parameters  == sde2.parameters
                             && sde1.periodicity == sde2.periodicity)

function Base.similar(sde::PSDE, q₀::AbstractArray, p₀::AbstractArray, ns::Int)
    similar(sde, sde.t₀, q₀, p₀, ns)
end

function Base.similar(sde::PSDE, q₀::AbstractArray, p₀::AbstractArray)
    similar(sde, sde.t₀, q₀, p₀)
end

function Base.similar(sde::PSDE, t₀::TT, q₀::AbstractArray{DT,N}, p₀::AbstractArray{DT,N}, ns::Int=sde.ns) where {DT <: Number, TT <: Number, N}
    @assert size(q₀) == size(p₀)
    @assert sde.d == size(q₀,1) == size(p₀,1)
    PSDE(sde.m, ns, sde.v, sde.f, sde.B, sde.G, t₀, q₀, p₀; parameters=sde.parameters, periodicity=sde.periodicity)
end

@inline Base.ndims(sde::PSDE) = sde.d

@inline periodicity(equation::PSDE) = equation.periodicity
