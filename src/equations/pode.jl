@doc raw"""
`PODE`: Partitioned Ordinary Differential Equation

Defines a partitioned initial value problem
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t), p(t)) , &
q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), p(t)) , &
p(t_{0}) &= p_{0} ,
\end{aligned}
```
with vector fields ``v`` and ``f``, initial conditions ``(q_{0}, p_{0})`` and the solution
``(q,p)`` taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}``.

### Parameters

* `DT <: Number`: data type
* `TT <: Real`: time step type
* `AT <: AbstractArray{DT}`: array type
* `vType <: Function`: type of `v`
* `fType <: Function`: type of `f`
* `hType <: OptionalFunction`: type of `h`
* `pType <: Union{NamedTuple,Nothing}`: parameters type

### Fields

* `d`: dimension of dynamical variables ``q`` and ``p`` as well as the vector fields ``v`` and ``f``
* `v`: function computing the vector field ``v``
* `f`: function computing the vector field ``f``
* `h`: function computing the Hamiltonian (optional)
* `t₀`: initial time
* `q₀`: initial condition for `q`
* `p₀`: initial condition for `p`
* `parameters`: either a `NamedTuple` containing the equations parameters or `nothing`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions

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

### Constructors

```julia
PODE(v, f, t₀, q₀, p₀; h=nothing, parameters=nothing, periodicity=zero(q₀[begin]))
PODE(v, f, q₀::StateVector, p₀::StateVector; kwargs...) = PODE(v, f, 0.0, q₀, p₀; kwargs...)
PODE(v, f, t₀, q₀::State, p₀::State; kwargs...) = PODE(v, f, t₀, [q₀], [p₀]; kwargs...)
PODE(v, f, q₀::State, p₀::State; kwargs...) = PODE(v, f, 0.0, q₀, p₀; kwargs...)
```

"""
struct PODE{dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
            vType <: Function, fType <: Function, hType <: OptionalFunction,
            pType <: Union{NamedTuple,Nothing}} <: AbstractEquationPODE{dType, tType}

    d::Int
    v::vType
    f::fType
    h::hType
    t₀::tType
    q₀::Vector{arrayType}
    p₀::Vector{arrayType}
    parameters::pType
    periodicity::arrayType

    function PODE(v::vType, f::fType, t₀::tType, q₀::Vector{arrayType}, p₀::Vector{arrayType};
                h::hType=nothing, parameters::pType=nothing,
                periodicity=zero(q₀[begin])) where {
                    dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
                    vType <: Function, fType <: Function,
                    hType <: OptionalFunction, pType <: Union{NamedTuple,Nothing}}

        d = length(q₀[begin])

        @assert length(q₀) == length(p₀)
        @assert all(length(q) == d for q in q₀)
        @assert all(length(p) == d for p in p₀)

        new{dType, tType, arrayType, vType, fType, hType, pType}(d, v, f, h, t₀, q₀, p₀, parameters, periodicity)
    end
end

PODE(v, f, q₀::StateVector, p₀::StateVector; kwargs...) = PODE(v, f, 0.0, q₀, p₀; kwargs...)
PODE(v, f, t₀, q₀::State, p₀::State; kwargs...) = PODE(v, f, t₀, [q₀], [p₀]; kwargs...)
PODE(v, f, q₀::State, p₀::State; kwargs...) = PODE(v, f, 0.0, q₀, p₀; kwargs...)

const PODEHT{HT,DT,TT,AT,VT,FT,PT} = PODE{DT,TT,AT,VT,FT,HT,PT} # type alias for dispatch on Hamiltonian type parameter
const PODEPT{PT,DT,TT,AT,VT,FT,HT} = PODE{DT,TT,AT,VT,FT,HT,PT} # type alias for dispatch on parameters type parameter


Base.hash(ode::PODE, h::UInt) = hash(ode.d, hash(ode.v, hash(ode.f, hash(ode.h,
        hash(ode.t₀, hash(ode.q₀, hash(ode.p₀, hash(ode.periodicity, hash(ode.parameters, h)))))))))

Base.:(==)(ode1::PODE, ode2::PODE) = (
                                ode1.d == ode2.d
                             && ode1.v == ode2.v
                             && ode1.f == ode2.f
                             && ode1.h == ode2.h
                             && ode1.t₀ == ode2.t₀
                             && ode1.q₀ == ode2.q₀
                             && ode1.p₀ == ode2.p₀
                             && ode1.parameters == ode2.parameters
                             && ode1.periodicity == ode2.periodicity)

Base.similar(equ::PODE, q₀, p₀; kwargs...) = similar(equ, equ.t₀, q₀, p₀; kwargs...)
Base.similar(equ::PODE, t₀::Real, q₀::State, p₀::State; kwargs...) = similar(equ, t₀, [q₀], [p₀]; kwargs...)

function Base.similar(equ::PODE, t₀::Real, q₀::StateVector, p₀::StateVector;
                      h=equ.h, parameters=equ.parameters, periodicity=equ.periodicity)
    @assert all([length(q) == ndims(equ) for q in q₀])
    @assert all([length(p) == ndims(equ) for p in p₀])
    PODE(equ.v, equ.f, t₀, q₀, p₀; h=h, parameters=parameters, periodicity=periodicity)
end

Base.ndims(equ::PODE) = equ.d
Base.axes(equ::PODE) = axes(equ.q₀[begin])
Common.nsamples(equ::PODE) = length(equ.q₀)
Common.periodicity(equation::PODE) = equation.periodicity

initial_conditions(equation::PODE) = (equation.t₀, equation.q₀, equation.p₀)

hashamiltonian(::PODEHT{<:Nothing}) = false
hashamiltonian(::PODEHT{<:Function}) = true

hasparameters(::PODEPT{<:Nothing}) = false
hasparameters(::PODEPT{<:NamedTuple}) = true


@define _create_pode_argument_views begin
    n = div(length(eachindex(x)), 2)
    q = @view x[eachindex(x)[  1:n ]]
    p = @view x[eachindex(x)[n+1:2n]]
    q̇ = @view ẋ[eachindex(ẋ)[  1:n ]]
    ṗ = @view ẋ[eachindex(ẋ)[n+1:2n]]
end


function Base.convert(::Type{ODE}, equ::PODE{DT,TT,AT}) where {DT, TT, AT <: AbstractVector}
    # concatenate initial conditions
    x₀ = [vcat(x...) for x in zip(equ.q₀, equ.p₀)]

    # extend periodicity
    periodicity = vcat(equ.periodicity, zero(equ.periodicity))

    if hasparameters(equ)
        v = (t, x, ẋ, params) -> begin
            @_create_pode_argument_views
            equ.v(t, q, p, q̇, params)
            equ.f(t, q, p, ṗ, params)
        end
    else
        v = (t, x, ẋ) -> begin
            @_create_pode_argument_views
            equ.v(t, q, p, q̇)
            equ.f(t, q, p, ṗ)
        end
    end

    ODE(v, equ.t₀, x₀; h=equ.h, parameters=equ.parameters, periodicity=periodicity)
end


_get_v(equ::PODE) = hasparameters(equ) ? (t,q,p,v) -> equ.v(t, q, p, v, equ.parameters) : equ.v
_get_f(equ::PODE) = hasparameters(equ) ? (t,q,p,f) -> equ.f(t, q, p, f, equ.parameters) : equ.f
_get_h(equ::PODE) = hasparameters(equ) ? (t,q,p) -> equ.h(t, q, p, equ.parameters) : equ.h
_get_v̄(equ::PODE) = _get_v(equ)
_get_f̄(equ::PODE) = _get_f(equ)


function get_function_tuple(equ::PODE)
    names = (:v,:f)
    equs  = (_get_v(equ), _get_f(equ))

    if hashamiltonian(equ)
        names = (names..., :h)
        equs  = (equs..., _get_h(equ))
    end

    NamedTuple{names}(equs)
end
