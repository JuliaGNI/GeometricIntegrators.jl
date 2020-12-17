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
* `q`: tuple of functions computing the solution
* `t₀`: initial time
* `q₀`: initial condition

The functions `v_i` providing the vector field must have the interface
```julia
    function v_i(t, q, v)
        v[1] = ...
        v[2] = ...
        ...
    end
```
and the functions `q_i` providing the solutions must have the interface
```julia
    function q_i(t, q₀, q₁, h)
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
struct SODE{dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
            vType <: Union{Tuple,Nothing}, qType <: Union{Tuple,Nothing},
            pType <: Union{NamedTuple,Nothing}} <: AbstractEquationODE{dType, tType}

    d::Int
    v::vType
    q::qType
    t₀::tType
    q₀::Vector{arrayType}
    parameters::pType
    periodicity::arrayType

    function SODE(v::vType, q::qType, t₀::tType, q₀::Vector{arrayType};
                 parameters::pType=nothing, periodicity=zero(q₀[begin])) where {
                        dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
                        vType <: Union{Tuple,Nothing}, qType <: Union{Tuple,Nothing},
                        pType <: Union{NamedTuple,Nothing}}

        d = length(q₀[begin])

        @assert all(length(q) == d for q in q₀)
        @assert v == nothing || all(typeof(V) <: OptionalFunction for V in v)
        @assert q == nothing || all(typeof(Q) <: OptionalFunction for Q in q)

        new{dType, tType, arrayType, vType, qType, pType}(d, v, q, t₀, q₀, parameters, periodicity)
    end
end

SODE(v, q::Union{Tuple,Nothing}, t₀::Real, q₀::State; kwargs...) = SODE(v, q, t₀, [q₀]; kwargs...)
SODE(v, q::Union{Tuple,Nothing}, q₀::StateVector; kwargs...) = SODE(v, q, 0.0, q₀; kwargs...)
SODE(v, q::Union{Tuple,Nothing}, q₀::State; kwargs...) = SODE(v, q, 0.0, q₀; kwargs...)

SODE(v, t₀::Real, q₀::StateVector; kwargs...) = SODE(v, nothing, t₀, q₀; kwargs...)
SODE(v, t₀::Real, q₀::State; kwargs...) = SODE(v, nothing, t₀, [q₀]; kwargs...)
SODE(v, q₀::StateVector; kwargs...) = SODE(v, nothing, 0.0, q₀; kwargs...)
SODE(v, q₀::State; kwargs...) = SODE(v, nothing, 0.0, q₀; kwargs...)

const SODEPT{PT,DT,TT,AT,VT,QT} = SODE{DT,TT,AT,VT,QT,PT} # type alias for dispatch on parameters type parameter
const SODEQT{QT,DT,TT,AT,VT,PT} = SODE{DT,TT,AT,VT,QT,PT} # type alias for dispatch on solution type parameter
const SODEVT{VT,DT,TT,AT,QT,PT} = SODE{DT,TT,AT,VT,QT,PT} # type alias for dispatch on vector field type parameter

Base.hash(ode::SODE, h::UInt) = hash(ode.d, hash(ode.v, hash(ode.t₀,
        hash(ode.q₀, hash(ode.periodicity, hash(ode.parameters, h))))))

Base.:(==)(ode1::SODE, ode2::SODE) = (
                                ode1.d == ode2.d
                             && ode1.v == ode2.v
                             && ode1.q == ode2.q
                             && ode1.t₀ == ode2.t₀
                             && ode1.q₀ == ode2.q₀
                             && ode1.parameters == ode2.parameters
                             && ode1.periodicity == ode2.periodicity)

Base.similar(equ::SODE, q₀; kwargs...) = similar(equ, equ.t₀, q₀; kwargs...)
Base.similar(equ::SODE, t₀::Real, q₀::State; kwargs...) = similar(equ, t₀, [q₀]; kwargs...)

function Base.similar(equ::SODE, t₀::Real, q₀::StateVector;
                      parameters=equ.parameters, periodicity=equ.periodicity)
    @assert all([length(q) == ndims(equ) for q in q₀])
    SODE(equ.v, equ.q, t₀, q₀; parameters=parameters, periodicity=periodicity)
end

Base.ndims(equ::SODE) = equ.d
Base.axes(equ::SODE) = axes(equ.q₀[begin])
Common.nsamples(equ::SODE) = length(equ.q₀)
Common.periodicity(equ::SODE) = equ.periodicity
initial_conditions(equ::SODE) = (equ.t₀, equ.q₀)

hasparameters(::SODEPT{<:Nothing}) = false
hasparameters(::SODEPT{<:NamedTuple}) = true

hassolution(::SODEQT{<:Nothing}) = false
hassolution(::SODEQT{<:Tuple}) = true # && all(typeof(Q) <: Functiong for Q in equ.q)

hassolution(::SODEQT{<:Nothing}, i) = false
hassolution(equ::SODEQT{<:Tuple}, i) = i ≤ length(equ.q) && typeof(equ.q[i]) <: Function

hasvectorfield(::SODEVT{<:Nothing}) = false
hasvectorfield(::SODEVT{<:Tuple}) = true # && all(typeof(V) <: Function for V in equ.v)

hasvectorfield(::SODEVT{<:Nothing}, i) = false
hasvectorfield(equ::SODEVT{<:Tuple}, i) = i ≤ length(equ.v) && typeof(equ.v[i]) <: Function

_get_v(equ::SODE) = hasparameters(equ) ? Tuple(typeof(V) <: Function ? (t,q,v) -> V(t, q, v, equ.parameters) : V for V in equ.v) : equ.v
_get_q(equ::SODE) = hasparameters(equ) ? Tuple(typeof(Q) <: Function ? (t,q̄,q,h) -> Q(t, q̄, q, h, equ.parameters) : Q for Q in equ.q) : equ.q

get_function_tuple(equ::SODE) = hasvectorfield(equ) ? _get_v(equ) : ()
get_solution_tuple(equ::SODE) = hassolution(equ) ? _get_q(equ) : ()
