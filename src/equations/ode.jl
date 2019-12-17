@doc raw"""
`ODE`: Ordinary Differential Equation

Defines an initial value problem
```math
\dot{q} (t) = v(t, q(t)) , \qquad q(t_{0}) = q_{0} ,
```
with vector field ``v``, initial condition ``q_{0}`` and the solution
``q`` taking values in ``\mathbb{R}^{d}``.

### Fields

* `d`: dimension of dynamical variable ``q`` and the vector field ``v``
* `n`: number of initial conditions
* `v`: function computing the vector field
* `t₀`: initial time
* `q₀`: initial condition

The function `v` providing the vector field must have the interface
```julia
    function v(t, q, v)
        v[1] = ...
        v[2] = ...
        ...
    end
```
where `t` is the current time, `q` is the current solution vector, and
`v` is the vector which holds the result of evaluating the vector field ``v``
on `t` and `q`.

"""
struct ODE{dType <: Number, tType <: Number,
           vType <: Function, pType <: Union{Tuple,Nothing}, N} <: AbstractEquationODE{dType, tType}

    d::Int
    n::Int
    v::vType
    t₀::tType
    q₀::Array{dType,N}
    parameters::pType
    periodicity::Vector{dType}

    function ODE(DT::DataType, N::Int, d::Int, n::Int, v::vType, t₀::tType, q₀::AbstractArray{dType};
                 parameters=nothing, periodicity=zeros(DT,d)) where {
                        dType <: Number, tType <: Number, vType <: Function}

        @assert d == size(q₀,1)
        @assert n == size(q₀,2)
        @assert ndims(q₀) == N ∈ (1,2)

        new{DT, tType, vType, typeof(parameters), N}(d, n, v, t₀,
                convert(Array{DT}, q₀), parameters, periodicity)
    end
end

function ODE(v, t₀, q₀::AbstractArray{DT}; kwargs...) where {DT}
    ODE(DT, ndims(q₀), size(q₀,1), size(q₀,2), v, t₀, q₀; kwargs...)
end

function ODE(v, q₀; kwargs...)
    ODE(v, zero(eltype(q₀)), q₀; kwargs...)
end

Base.hash(ode::ODE, h::UInt) = hash(ode.d, hash(ode.n, hash(ode.v, hash(ode.t₀,
        hash(ode.q₀, hash(ode.periodicity, hash(ode.parameters, h)))))))

Base.:(==)(ode1::ODE, ode2::ODE) = (
                                ode1.d == ode2.d
                             && ode1.n == ode2.n
                             && ode1.v == ode2.v
                             && ode1.t₀ == ode2.t₀
                             && ode1.q₀ == ode2.q₀
                             && ode1.parameters == ode2.parameters
                             && ode1.periodicity == ode2.periodicity)

function Base.similar(ode::ODE, q₀; kwargs...)
    similar(ode, ode.t₀, q₀; kwargs...)
end

function Base.similar(ode::ODE, t₀::TT, q₀::AbstractArray{DT};
                      parameters=ode.parameters, periodicity=ode.periodicity) where {DT  <: Number, TT <: Number}
    @assert ode.d == size(q₀,1)
    ODE(ode.v, t₀, q₀; parameters=parameters, periodicity=periodicity)
end

@inline Base.ndims(ode::ODE) = ode.d

@inline periodicity(equation::ODE) = equation.periodicity
