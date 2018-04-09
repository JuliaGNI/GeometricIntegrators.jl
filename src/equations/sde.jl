
"""
`SDE`: Stratonovich Stochastic Differential Equation

Defines a stochastic differential initial value problem
```math
\\begin{align*}
\\dq (t) &= v(t, q(t)) \, dt + B(t, q(t)) \circ dW , & q(t_{0}) &= q_{0} ,
\\end{align*}
```
with drift vector field ``v``, diffusion matrix ``B``,
initial conditions ``q_{0}``, the dynamical variable ``q``
taking values in ``\\mathbb{R}^{d}``, and the m-dimensional Wiener process W

### Fields

* `d`:  dimension of dynamical variable ``q`` and the vector field ``v``
* `m`:  dimension of the Wiener process
* `n`:  number of initial conditions
* `ns`: number of sample paths
* `v`:  function computing the deterministic vector field
* `B`:  function computing the d x m diffusion matrix
* `t₀`: initial time
* `q₀`: initial condition for dynamical variable ``q`` (may be a random variable itself)


The functions `v` and `B`, providing the drift vector field and diffusion matrix, take three arguments,
`v(t, q, v)` and `B(t, q, B)`, where `t` is the current time, `q` is the
current solution vector, and `v` and `B` are the variables which hold the result
of evaluating the vector field ``v`` and the matrix ``B`` on `t` and `q`.

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
struct SDE{dType <: Number, tType <: Number, vType <: Function, BType <: Function, N} <: Equation{dType, tType}
    d::Int
    m::Int
    n::Int
    ns::Int
    v::vType
    B::BType
    t₀::tType
    q₀::Array{dType, N}           #Initial condition: N=1 - single deterministic, N=2 - single random or multiple deterministic, N=3 - multiple deterministic
    periodicity::Vector{dType}

    function SDE{dType,tType,vType,BType,N}(d, m, n, ns, v, B, t₀, q₀; periodicity=[]) where {dType <: Number, tType <: Number, vType <: Function, BType <: Function, N}

        @assert dType == eltype(q₀)
        @assert tType == typeof(t₀)
        @assert ndims(q₀) == N ∈ (1,2,3)
        @assert d == size(q₀,1)

        if ns==1 && n==1
            @assert N == 1
        elseif ns>1 && n==1
            @assert N == 2
            @assert ns == size(q₀,2)
        elseif ns==1 && n>1
            @assert N == 2
            @assert n == size(q₀,2)
        elseif ns>1 && n>1
            @assert N == 3
            @assert ns == size(q₀,2)
            @assert n == size(q₀,3)
        end

        if !(length(periodicity) == d)
            periodicity = zeros(dType, d)
        end

        new(d, m, n, ns, v, B, t₀, q₀, periodicity)
    end
end


function SDE(m::Int, v::VT, B::BT, t₀::TT, q₀::DenseArray{DT,1}; periodicity=[]) where {DT,TT,VT,BT}
    SDE{DT, TT, VT, BT, 1}(size(q₀, 1), m, 1, 1, v, B, t₀, q₀, periodicity=periodicity)
end


# A 2-dimensional matrix q0 can represent a single random initial condition with ns>1 and n=1,
# or a set of deterministic initial conditions with ns=1 and n>1.
# The argument IC specifies whether there are multiple deterministic initial conditions
# (IC=true, so n>1) or a single random one (default IC=false, so n=1)
function SDE(m::Int, v::VT, B::BT, t₀::TT, q₀::DenseArray{DT,2}; periodicity=[], IC=false) where {DT,TT,VT,BT}
    if IC==true
        SDE{DT, TT, VT, BT, 2}(size(q₀, 1), m, size(q₀, 2), 1, v, B, t₀, q₀, periodicity=periodicity)
    else
        SDE{DT, TT, VT, BT, 2}(size(q₀, 1), m, 1, size(q₀, 2), v, B, t₀, q₀, periodicity=periodicity)
    end
end


function SDE(m::Int, v::VT, B::BT, t₀::TT, q₀::DenseArray{DT,3}; periodicity=[]) where {DT,TT,VT,BT}
    SDE{DT, TT, VT, BT, 3}(size(q₀, 1), m, size(q₀,3), size(q₀,2), v, B, t₀, q₀, periodicity=periodicity)
end


function SDE(m::Int, v::VT, B::BT, q₀::DenseArray{DT,1}; periodicity=[]) where {DT,VT,BT}
    SDE(m, v, B, zero(DT), q₀, periodicity=periodicity)
end


function SDE(m::Int, v::VT, B::BT, q₀::DenseArray{DT,2}; periodicity=[], IC=false) where {DT,VT,BT}
    SDE(m, v, B, zero(DT), q₀, periodicity=periodicity, IC=IC)
end


function SDE(m::Int, v::VT, B::BT, q₀::DenseArray{DT,3}; periodicity=[]) where {DT,VT,BT}
    SDE(m, v, B, zero(DT), q₀, periodicity=periodicity)
end

Base.hash(sde::SDE, h::UInt) = hash(sde.d, hash(sde.m, hash(sde.n, hash(sde.ns, hash(sde.v, hash(sde.B, hash(sde.t₀, hash(sde.q₀, hash(sde.periodicity, h)))))))))

Base.:(==)(sde1::SDE, sde2::SDE) = (
                                sde1.d == sde2.d
                             && sde1.m == sde2.m
                             && sde1.n == sde2.n
                             && sde1.ns == sde2.ns
                             && sde1.v == sde2.v
                             && sde1.B == sde2.B
                             && sde1.t₀ == sde2.t₀
                             && sde1.q₀ == sde2.q₀
                             && sde1.periodicity == sde2.periodicity)

function Base.similar(sde::SDE{DT,TT,VT,UT}, q₀::DenseArray{DT,1}) where {DT, TT, VT, UT}
    similar(sde, sde.t₀, q₀)
end

function Base.similar(sde::SDE{DT,TT,VT,UT}, q₀::DenseArray{DT,2}; IC=false) where {DT, TT, VT, UT}
    similar(sde, sde.t₀, q₀, IC=IC)
end

function Base.similar(sde::SDE{DT,TT,VT,UT}, q₀::DenseArray{DT,3}) where {DT, TT, VT, UT}
    similar(sde, sde.t₀, q₀)
end

function Base.similar(sde::SDE{DT,TT,VT,BT}, t₀::TT, q₀::DenseArray{DT,1}) where {DT, TT, VT, BT}
    @assert sde.d == size(q₀,1)
    SDE(sde.m, sde.v, sde.B, t₀, q₀, periodicity=sde.periodicity)
end

function Base.similar(sde::SDE{DT,TT,VT,BT}, t₀::TT, q₀::DenseArray{DT,2}; IC=false) where {DT, TT, VT, BT}
    @assert sde.d == size(q₀,1)
    SDE(sde.m, sde.v, sde.B, t₀, q₀, periodicity=sde.periodicity, IC=IC)
end

function Base.similar(sde::SDE{DT,TT,VT,BT}, t₀::TT, q₀::DenseArray{DT,3}) where {DT, TT, VT, BT}
    @assert sde.d == size(q₀,1)
    SDE(sde.m, sde.v, sde.B, t₀, q₀, periodicity=sde.periodicity)
end
