@doc raw"""
`PODE`: Partitioned Ordinary Differential Equation

Defines a partitioned initial value problem
```math
\begin{align*}
\dot{q} (t) &= v(t, q(t), p(t)) , &
q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), p(t)) , &
p(t_{0}) &= p_{0} ,
\end{align*}
```
with vector fields ``v`` and ``f``, initial conditions ``(q_{0}, p_{0})`` and the solution
``(q,p)`` taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}``.

### Fields

* `d`: dimension of dynamical variables ``q`` and ``p`` as well as the vector fields ``v`` and ``f``
* `v`: function computing the vector field ``v``
* `f`: function computing the vector field ``f``
* `t₀`: initial time
* `q₀`: initial condition for `q`
* `p₀`: initial condition for `p`

The functions `v` and `f` must have the interface
```julia
    function v(t, q, p, v)
        v[1] = ...
        v[2] = ...
        ...
    end
```
and
```julia
    function f(t, q, p, f)
        f[1] = ...
        f[2] = ...
        ...
    end
```
where `t` is the current time, `q` and `p` are the current solution vectors
and `v` and `f` are the vectors which hold the result of evaluating the
vector fields ``v`` and ``f`` on `t`, `q` and `p`.
"""
struct PODE{dType <: Number, tType <: Number,
            vType <: Function, fType <: Function,
            pType <: Union{Tuple,Nothing}, N} <: AbstractEquationPODE{dType, tType}

    d::Int
    n::Int
    v::vType
    f::fType
    t₀::tType
    q₀::Array{dType, N}
    p₀::Array{dType, N}
    parameters::pType
    periodicity::Vector{dType}

    function PODE(DT::DataType, N::Int, d::Int, n::Int, v::vType, f::fType,
                  t₀::tType, q₀::AbstractArray{dType}, p₀::AbstractArray{dType};
                  parameters=nothing, periodicity=zeros(DT,d)) where {
                        dType <: Number, tType <: Number, vType <: Function, fType <: Function}

        @assert d == size(q₀,1) == size(p₀,1)
        @assert n == size(q₀,2) == size(p₀,2)
        @assert ndims(q₀) == ndims(p₀) == N ∈ (1,2)

        new{DT, tType, vType, fType, typeof(parameters), N}(d, n, v, f, t₀,
                convert(Array{DT}, q₀), convert(Array{DT}, p₀),
                parameters, periodicity)
    end
end

function PODE(v, f, t₀, q₀::AbstractArray{DT}, p₀::AbstractArray{DT}; kwargs...) where {DT}
    PODE(DT, ndims(q₀), size(q₀,1), size(q₀,2), v, f, t₀, q₀, p₀; kwargs...)
end

function PODE(v, f, q₀, p₀; kwargs...)
    PODE(v, f, zero(eltype(q₀)), q₀, p₀; kwargs...)
end

Base.hash(ode::PODE, h::UInt) = hash(ode.d, hash(ode.n, hash(ode.v, hash(ode.f,
        hash(ode.t₀, hash(ode.q₀, hash(ode.p₀, hash(ode.periodicity, hash(ode.parameters, h)))))))))

Base.:(==)(ode1::PODE{DT1,TT1,VT1,FT1}, ode2::PODE{DT2,TT2,VT2,FT2}) where {DT1, DT2, TT1, TT2, VT1, VT2, FT1, FT2} = (
                                ode1.d == ode2.d
                             && ode1.n == ode2.n
                             && ode1.v == ode2.v
                             && ode1.f == ode2.f
                             && ode1.t₀ == ode2.t₀
                             && ode1.q₀ == ode2.q₀
                             && ode1.p₀ == ode2.p₀
                             && ode1.parameters == ode2.parameters
                             && ode1.periodicity == ode2.periodicity)

function Base.similar(ode::PODE, q₀, p₀; kwargs...)
    similar(ode, ode.t₀, q₀, p₀; kwargs...)
end

function Base.similar(ode::PODE, t₀::TT, q₀::AbstractArray{DT}, p₀::AbstractArray{DT};
                      parameters=ode.parameters, periodicity=ode.periodicity) where {DT  <: Number, TT <: Number}
    @assert ode.d == size(q₀,1) == size(p₀,1)
    PODE(ode.v, ode.f, t₀, q₀, p₀; parameters=parameters, periodicity=periodicity)
end

@inline Base.ndims(ode::PODE) = ode.d

@inline CommonFunctions.periodicity(equation::PODE) = equation.periodicity
