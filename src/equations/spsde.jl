@doc raw"""
`SPSDE`: Stratonovich Split Partitioned Stochastic Differential Equation

Defines a partitioned stochastic differential initial value problem
```math
\begin{align*}
\dq (t) &=   v(t, q(t)) \, dt + B(t, q(t)) \circ dW , & q(t_{0}) &= q_{0} ,
\dp (t) &= [ f1(t, q(t)) + f2(t, q(t)) ] \, dt + [ G1(t, q(t)) + G2(t, q(t)) ] \circ dW , & p(t_{0}) &= p_{0}
\end{align*}
```
with the drift vector fields ``v`` and ``fi``, diffusion matrices ``B`` and ``Gi``,
initial conditions ``q_{0}`` and ``p_{0}``, the dynamical variables ``(q,p)`` taking
values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}``, and the m-dimensional Wiener process W

### Fields

* `d`:  dimension of dynamical variable ``q`` and the vector fields ``vi``
* `m`:  dimension of the Wiener process
* `n`:  number of initial conditions
* `ns`: number of sample paths
* `v` :  function computing the drift vector field for the position variable q
* `f1`:  function computing the drift vector field for the momentum variable p
* `f2`:  function computing the drift vector field for the momentum variable p
* `B` :  function computing the d x m diffusion matrix for the position variable q
* `G1`:  function computing the d x m diffusion matrix for the momentum variable p
* `G2`:  function computing the d x m diffusion matrix for the momentum variable p
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
struct SPSDE{dType <: Number, tType <: Number, vType <: Function, f1Type <: Function, f2Type <: Function, BType <: Function, G1Type <: Function, G2Type <: Function, N} <: Equation{dType, tType}
    d::Int
    m::Int
    n::Int
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
    periodicity::Vector{dType}

    function SPSDE{dType,tType,vType,f1Type,f2Type,BType,G1Type,G2Type,N}(d, m, n, ns, v, f1, f2, B, G1, G2, t₀, q₀, p₀; periodicity=[]) where {dType <: Number, tType <: Number, vType <: Function, f1Type <: Function, f2Type <: Function, BType <: Function, G1Type <: Function, G2Type <: Function, N}

        @assert dType == eltype(q₀) == eltype(p₀)
        @assert tType == typeof(t₀)
        @assert ndims(q₀) == ndims(p₀) == N
        @assert d == size(q₀,1) == size(p₀,1)

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

        if !(length(periodicity) == d)
            periodicity = zeros(dType, d)
        end

        new(d, m, n, ns, v, f1, f2, B, G1, G2, t₀, q₀, p₀, periodicity)
    end
end


function SPSDE(m::Int, ns::Int, v::VT, f1::F1T, f2::F2T, B::BT, G1::G1T, G2::G2T, t₀::TT, q₀::DenseArray{DT,1}, p₀::DenseArray{DT,1}; periodicity=[]) where {DT,TT,VT,F1T,F2T,BT,G1T,G2T}
    # A 1D array q₀ contains a single deterministic initial condition, so n=1, but we still need to specify
    # the number of sample paths ns
    @assert size(q₀) == size(p₀)
    SPSDE{DT, TT, VT, F1T, F2T, BT, G1T, G2T, 1}(size(q₀, 1), m, 1, ns, v, f1, f2, B, G1, G2, t₀, q₀, p₀, periodicity=periodicity)
end


# A 2-dimensional matrix q0 can represent a single random initial condition with ns>1 and n=1,
# or a set of deterministic initial conditions with n>1 (for which we can have both ns=1 and ns>1)
# The function below assumes q0 to represent a single random initial condition (n=1, ns=size(q₀, 2))
function SPSDE(m::Int, v::VT, f1::F1T, f2::F2T, B::BT, G1::G1T, G2::G2T, t₀::TT, q₀::DenseArray{DT,2}, p₀::DenseArray{DT,2}; periodicity=[]) where {DT,TT,VT,F1T,F2T,BT,G1T,G2T}
    SPSDE{DT, TT, VT, F1T, F2T, BT, G1T, G2T, 2}(size(q₀, 1), m, 1, size(q₀, 2), v, f1, f2, B, G1, G2, t₀, q₀, p₀, periodicity=periodicity)
end

# On the other hand, the function below assumes q₀ represents multiple deterministic initial conditions
# (n=size(q₀, 2)), but these initial conditions may be run an arbitrary number ns of sample paths, so ns has to be explicitly specified
function SPSDE(m::Int, ns::Int, v::VT, f1::F1T, f2::F2T, B::BT, G1::G1T, G2::G2T, t₀::TT, q₀::DenseArray{DT,2}, p₀::DenseArray{DT,2}; periodicity=[]) where {DT,TT,VT,F1T,F2T,BT,G1T,G2T}
    SPSDE{DT, TT, VT, F1T, F2T, BT, G1T, G2T, 2}(size(q₀, 1), m, size(q₀, 2), ns, v, f1, f2, B, G1, G2, t₀, q₀, p₀, periodicity=periodicity)
end


function SPSDE(m::Int, v::VT, f1::F1T, f2::F2T, B::BT, G1::G1T, G2::G2T, t₀::TT, q₀::DenseArray{DT,3}, p₀::DenseArray{DT,3}; periodicity=[]) where {DT,TT,VT,F1T,F2T,BT,G1T,G2T}
    @assert size(q₀) == size(p₀)
    SPSDE{DT, TT, VT, F1T, F2T, BT, G1T, G2T, 3}(size(q₀, 1), m, size(q₀,3), size(q₀,2), v, f1, f2, B, G1, G2, t₀, q₀, p₀, periodicity=periodicity)
end


function SPSDE(m::Int, ns::Int, v::VT, f1::F1T, f2::F2T, B::BT, G1::G1T, G2::G2T, q₀::DenseArray{DT,1}, p₀::DenseArray{DT,1}; periodicity=[]) where {DT,VT,F1T,F2T,BT,G1T,G2T}
    SPSDE(m, ns, v, f1, f2, B, G1, G2, zero(DT), q₀, p₀, periodicity=periodicity)
end


# Assumes q0 represents a single random initial condition (n=1, ns=size(q₀, 2))
function SPSDE(m::Int, v::VT, f1::F1T, f2::F2T, B::BT, G1::G1T, G2::G2T, q₀::DenseArray{DT,2}, p₀::DenseArray{DT,2}; periodicity=[]) where {DT,VT,F1T,F2T,BT,G1T,G2T}
    SPSDE(m, v, f1, f2, B, G1, G2, zero(DT), q₀, p₀, periodicity=periodicity)
end


# Assumes q₀ represents multiple deterministic initial conditions (n=size(q₀, 2))
function SPSDE(m::Int, ns::Int, v::VT, f1::F1T, f2::F2T, B::BT, G1::G1T, G2::G2T, q₀::DenseArray{DT,2}, p₀::DenseArray{DT,2}; periodicity=[]) where {DT,VT,F1T,F2T,BT,G1T,G2T}
    SPSDE(m, ns, v, f1, f2, B, G1, G2, zero(DT), q₀, p₀, periodicity=periodicity)
end


function SPSDE(m::Int, v::VT, f1::F1T, f2::F2T, B::BT, G1::G1T, G2::G2T, q₀::DenseArray{DT,3}, p₀::DenseArray{DT,3}; periodicity=[]) where {DT,VT,F1T,F2T,BT,G1T,G2T}
    SPSDE(m, v, f1, f2, B, G1, G2, zero(DT), q₀, p₀, periodicity=periodicity)
end

Base.hash(sde::SPSDE, h::UInt) = hash(sde.d, hash(sde.m, hash(sde.n, hash(sde.ns, hash(sde.v, hash(sde.f1, hash(sde.f2, hash(sde.B, hash(sde.G1, hash(sde.G2, hash(sde.t₀, hash(sde.q₀, hash(sde.p₀, hash(sde.periodicity, h))))))))))))))

Base.:(==)(sde1::SPSDE, sde2::SPSDE) = (
                                sde1.d == sde2.d
                             && sde1.m == sde2.m
                             && sde1.n == sde2.n
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
                             && sde1.periodicity == sde2.periodicity)

function Base.similar(sde::SPSDE{DT,TT,VT,F1T,F2T,BT,G1T,G2T}, q₀::DenseArray{DT,1}, p₀::DenseArray{DT,1}, ns::Int) where {DT, TT, VT, F1T, F2T, BT, G1T, G2T}
    similar(sde, sde.t₀, q₀, p₀, ns)
end

# Assumes q0 represents a single random initial condition (n=1, ns=size(q₀, 2))
function Base.similar(sde::SPSDE{DT,TT,VT,F1T,F2T,BT,G1T,G2T}, q₀::DenseArray{DT,2}, p₀::DenseArray{DT,2}) where {DT, TT, VT, F1T, F2T, BT, G1T, G2T}
    similar(sde, sde.t₀, q₀, p₀)
end

# Assumes q₀ represents multiple deterministic initial conditions (n=size(q₀, 2))
function Base.similar(sde::SPSDE{DT,TT,VT,F1T,F2T,BT,G1T,G2T}, q₀::DenseArray{DT,2}, p₀::DenseArray{DT,2}, ns::Int) where {DT, TT, VT, F1T, F2T, BT, G1T, G2T}
    similar(sde, sde.t₀, q₀, p₀, ns)
end

function Base.similar(sde::SPSDE{DT,TT,VT,F1T,F2T,BT,G1T,G2T}, q₀::DenseArray{DT,3}, p₀::DenseArray{DT,3}) where {DT, TT, VT, F1T, F2T, BT, G1T, G2T}
    similar(sde, sde.t₀, q₀, p₀)
end

function Base.similar(sde::SPSDE{DT,TT,VT,F1T,F2T,BT,G1T,G2T}, t₀::TT, q₀::DenseArray{DT,1}, p₀::DenseArray{DT,1}, ns::Int) where {DT, TT, VT, F1T, F2T, BT, G1T, G2T}
    @assert size(q₀) == size(p₀)
    @assert sde.d == size(q₀,1)
    SPSDE(sde.m, ns, sde.v, sde.f1, sde.f2, sde.B, sde.G1, sde.G2, t₀, q₀, p₀, periodicity=sde.periodicity)
end

# Assumes q0 represents a single random initial condition (n=1, ns=size(q₀, 2))
function Base.similar(sde::SPSDE{DT,TT,VT,F1T,F2T,BT,G1T,G2T}, t₀::TT, q₀::DenseArray{DT,2}, p₀::DenseArray{DT,2}) where {DT, TT, VT, F1T, F2T, BT, G1T, G2T}
    @assert size(q₀) == size(p₀)
    @assert sde.d == size(q₀,1)
    SPSDE(sde.m, sde.v, sde.f1, sde.f2, sde.B, sde.G1, sde.G2, t₀, q₀, p₀, periodicity=sde.periodicity)
end

# Assumes q₀ represents multiple deterministic initial conditions (n=size(q₀, 2))
function Base.similar(sde::SPSDE{DT,TT,VT,F1T,F2T,BT,G1T,G2T}, t₀::TT, q₀::DenseArray{DT,2}, p₀::DenseArray{DT,2}, ns::Int) where {DT, TT, VT, F1T, F2T, BT, G1T, G2T}
    @assert size(q₀) == size(p₀)
    @assert sde.d == size(q₀,1)
    SPSDE(sde.m, ns, sde.v, sde.f1, sde.f2, sde.B, sde.G1, sde.G2, t₀, q₀, p₀, periodicity=sde.periodicity)
end

function Base.similar(sde::SPSDE{DT,TT,VT,F1T,F2T,BT,G1T,G2T}, t₀::TT, q₀::DenseArray{DT,3}, p₀::DenseArray{DT,3}) where {DT, TT, VT, F1T, F2T, BT, G1T, G2T}
    @assert size(q₀) == size(p₀)
    @assert sde.d == size(q₀,1)
    PSDE(sde.m, sde.v, sde.f1, sde.f2, sde.B, sde.G1, sde.G2, t₀, q₀, p₀, periodicity=sde.periodicity)
end

@inline Base.ndims(sde::SPSDE) = sde.d

@inline periodicity(equation::SPSDE) = equation.periodicity
