@doc raw"""
`SODE`: Split Ordinary Differential Equation

Defines an initial value problem
```math
\dot{q} (t) = v(t, q(t)) , \qquad q(t_{0}) = q_{0} ,
```
with vector field ``v``, initial condition ``q_{0}`` and the solution
``q`` taking values in ``\mathbb{R}^{d}``. Here, the vector field ``v``
is given as a sum of vector fields
```math
v (t) = v_1 (t) + ... + v_r (t) .
```

### Fields

* `d`: dimension of dynamical variable ``q`` and the vector field ``v``
* `v`: tuple of functions computing the vector field
* `t₀`: initial time
* `q₀`: initial condition

The functions `v_i` providing the vector field must have the interface
```julia
    function v_i(t, q₀, q₁, h)
        q₁[1] = q₀[1] + ...
        q₁[2] = q₀[2] + ...
        ...
    end
```
where `t` is the current time, `q₀` is the current solution vector, `q₁` is the
new solution vector which holds the result of computing one substep with the
vector field ``v_i`` on `t` and `q₀`, and `h` is the (sub-)timestep to compute
the update for.

The fact that the function `v` returns the solution and not just the vector
field for each substep increases the flexibility for the use of splitting
methods, e.g., it allows to use another integrator for solving substeps.

"""
struct SODE{dType <: Number, tType <: Number,
            vType <: Tuple, pType <: Union{Tuple,Nothing}, N} <: AbstractEquationODE{dType, tType}

    d::Int
    n::Int
    v::vType
    t₀::tType
    q₀::Array{dType,N}
    parameters::pType
    periodicity::Vector{dType}

    function SODE(DT::DataType, N::Int, d::Int, n::Int, v::vType, t₀::tType, q₀::AbstractArray{dType};
                 parameters=nothing, periodicity=zeros(DT,d)) where {
                        dType <: Number, tType <: Number, vType <: Tuple}

        @assert d == size(q₀,1)
        @assert n == size(q₀,2)
        @assert ndims(q₀) == N ∈ (1,2)

        new{DT, tType, vType, typeof(parameters), N}(d, n, v, t₀,
                convert(Array{DT}, q₀), parameters, periodicity)
    end
end

function SODE(v, t₀, q₀::AbstractArray{DT}; kwargs...) where {DT}
    SODE(DT, ndims(q₀), size(q₀,1), size(q₀,2), v, t₀, q₀; kwargs...)
end

function SODE(v, q₀; kwargs...)
    SODE(v, zero(eltype(q₀)), q₀; kwargs...)
end

Base.hash(ode::SODE, h::UInt) = hash(ode.d, hash(ode.n, hash(ode.v, hash(ode.t₀,
        hash(ode.q₀, hash(ode.periodicity, hash(ode.parameters, h)))))))

Base.:(==)(ode1::SODE, ode2::SODE) = (
                                ode1.d == ode2.d
                             && ode1.n == ode2.n
                             && ode1.v == ode2.v
                             && ode1.t₀ == ode2.t₀
                             && ode1.q₀ == ode2.q₀
                             && ode1.parameters == ode2.parameters
                             && ode1.periodicity == ode2.periodicity)

function Base.similar(ode::SODE, q₀; kwargs...)
    similar(ode, ode.t₀, q₀; kwargs...)
end

function Base.similar(ode::SODE, t₀::TT, q₀::AbstractArray{DT};
                      parameters=ode.parameters, periodicity=ode.periodicity) where {DT  <: Number, TT <: Number}
    @assert ode.d == size(q₀,1)
    SODE(ode.v, t₀, q₀; parameters=parameters, periodicity=periodicity)
end

@inline Base.ndims(ode::SODE) = ode.d

@inline CommonFunctions.periodicity(equation::SODE) = equation.periodicity
