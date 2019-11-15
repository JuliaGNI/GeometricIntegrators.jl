@doc raw"""
`SDE`: Stratonovich Stochastic Differential Equation

Defines a stochastic differential initial value problem
```math
\begin{align*}
\dq (t) &= v(t, q(t)) \, dt + B(t, q(t)) \circ dW , & q(t_{0}) &= q_{0} ,
\end{align*}
```
with drift vector field ``v``, diffusion matrix ``B``,
initial conditions ``q_{0}``, the dynamical variable ``q``
taking values in ``\mathbb{R}^{d}``, and the m-dimensional Wiener process W

### Fields

* `d`:  dimension of dynamical variable ``q`` and the vector field ``v``
* `m`:  dimension of the Wiener process
* `n`:  number of initial conditions
* `ns`: number of sample paths
* `v`:  function computing the deterministic vector field
* `B`:  function computing the d x m diffusion matrix
* `t₀`: initial time
* `q₀`: initial condition for dynamical variable ``q`` (may be a random variable itself)


The functions `v` and `B`, providing the drift vector field and diffusion matrix,
`v(t, q, v)` and `B(t, q, B; col=0)`, where `t` is the current time, `q` is the
current solution vector, and `v` and `B` are the variables which hold the result
of evaluating the vector field ``v`` and the matrix ``B`` on `t` and `q` (if col==0),
or the column col of the matrix B (if col>0).

### Example

```julia
    function v(λ, t, q, v)
        v[1] = λ*q[1]
        v[2] = λ*q[2]
    end

    function B(μ, t, q, B; col=0)
        if col==0 #whole matrix
            B[1,1] = μ*q[1]
            B[2,1] = μ*q[2]
        elseif col==1
            #just first column
        end
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
struct SDE{dType <: Number, tType <: Number, vType <: Function, BType <: Function,
           pType <: Union{Tuple,Nothing}, N} <: Equation{dType, tType}

    d::Int
    m::Int
    n::Int
    ns::Int
    v::vType
    B::BType
    t₀::tType
    q₀::Array{dType, N}           #Initial condition: N=1 - single deterministic, N=2 - single random or multiple deterministic, N=3 - multiple random
    parameters::pType
    periodicity::Vector{dType}

    function SDE(m, n, ns, v::vType, B::BType, t₀::tType, q₀::DenseArray{dType};
                 parameters=nothing, periodicity=zeros(dType,size(q₀,1))) where {
                        dType <: Number, tType <: Number, vType <: Function, BType <: Function}

        N = ndims(q₀)
        d = size(q₀,1)

        if ns==1 && n==1
            # single sample path and single initial condition, therefore N=1
            @assert N == 1
        elseif ns>1 && n==1
            # single initial condition, but may be random (N=2) or deterministic (N=1)
            @assert N ∈ (1,2)
            if N==2
                @assert ns == size(q₀,2)
            end
        elseif ns==1 && n>1
            # multiple  deterministic initial conditions, so N=2
            @assert N == 2
            @assert n == size(q₀,2)
        elseif ns>1 && n>1
            # either multiple random initial conditions (N=3),
            # or multiple deterministic initial conditions (N=2)
            @assert N ∈ (2,3)

            if N==2
                @assert n == size(q₀,2)
            else
                @assert ns == size(q₀,2)
                @assert n == size(q₀,3)
            end
        end

        new{dType,tType,vType,BType,typeof(parameters),N}(d, m, n, ns, v, B, t₀, q₀, parameters, periodicity)
    end
end


function SDE(m::Int, ns::Int, v::Function, B::Function, t₀::Number, q₀::DenseArray{DT,1}; kwargs...) where {DT <: Number}
    # A 1D array q₀ contains a single deterministic initial condition, so n=1, but we still need to specify
    # the number of sample paths ns
    SDE(m, 1, ns, v, B, t₀, q₀; kwargs...)
end

# A 2-dimensional matrix q0 can represent a single random initial condition with ns>1 and n=1,
# or a set of deterministic initial conditions with n>1 (for which we can have both ns=1 and ns>1)
# The function below assumes q0 to represent a single random initial condition (n=1, ns=size(q₀, 2))
function SDE(m::Int, v::Function, B::Function, t₀::Number, q₀::DenseArray{DT,2}; kwargs...) where {DT <: Number}
    SDE(m, 1, size(q₀,2), v, B, t₀, q₀; kwargs...)
end

# On the other hand, the function below assumes q₀ represents multiple deterministic initial conditions
# (n=size(q₀, 2)), but these initial conditions may be run an arbitrary number ns of sample paths, so ns has to be explicitly specified
function SDE(m::Int, ns::Int, v::Function, B::Function, t₀::Number, q₀::DenseArray{DT,2}; kwargs...) where {DT <: Number}
    SDE(m, size(q₀,2), ns, v, B, t₀, q₀; kwargs...)
end

function SDE(m::Int, v::Function, B::Function, t₀::Number, q₀::DenseArray{DT,3}; kwargs...) where {DT <: Number}
    # A 3D array q₀ contains multiple random initial condition, so n=size(q₀,3) and ns=size(q₀,2)
    SDE(m, size(q₀,3), size(q₀,2), v, B, t₀, q₀; kwargs...)
end


function SDE(m::Int, ns::Int, v::Function, B::Function, q₀::DenseArray{DT}; kwargs...) where {DT}
    SDE(m, ns, v, B, zero(DT), q₀; kwargs...)
end

function SDE(m::Int, v::Function, B::Function, q₀::DenseArray{DT}; kwargs...) where {DT}
    SDE(m, v, B, zero(DT), q₀; kwargs...)
end


Base.hash(sde::SDE, h::UInt) = hash(sde.d, hash(sde.m, hash(sde.n, hash(sde.ns,
                               hash(sde.v, hash(sde.B, hash(sde.t₀, hash(sde.q₀,
                               hash(sde.parameters, hash(sde.periodicity, h))))))))))

Base.:(==)(sde1::SDE, sde2::SDE) = (
                                sde1.d == sde2.d
                             && sde1.m == sde2.m
                             && sde1.n == sde2.n
                             && sde1.ns == sde2.ns
                             && sde1.v == sde2.v
                             && sde1.B == sde2.B
                             && sde1.t₀ == sde2.t₀
                             && sde1.q₀ == sde2.q₀
                             && sde1.parameters == sde2.parameters
                             && sde1.periodicity == sde2.periodicity)

function Base.similar(sde::SDE, q₀::DenseArray)
    similar(sde, sde.t₀, q₀)
end

function Base.similar(sde::SDE, t₀::TT, q₀::DenseArray{DT}) where {DT <: Number, TT <: Number}
    @assert sde.d == size(q₀,1)
    SDE(sde.m, sde.ns, sde.v, sde.B, t₀, q₀, periodicity=sde.periodicity)
end

function Base.similar(sde::SDE, q₀::DenseArray, ns::Int)
    similar(sde, sde.t₀, q₀, ns)
end

function Base.similar(sde::SDE, t₀::TT, q₀::DenseArray{DT}, ns::Int) where {DT <: Number, TT <: Number}
    @assert sde.d == size(q₀,1)
    SDE(sde.m, ns, sde.v, sde.B, t₀, q₀, periodicity=sde.periodicity)
end

Base.ndims(sde::SDE) = sde.d
