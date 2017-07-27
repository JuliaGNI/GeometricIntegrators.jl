"""
`SODE`: Split Ordinary Differential Equation

Defines an initial value problem
```math
\\dot{q} (t) = v(t, q(t)) , \\qquad q(t_{0}) = q_{0} ,
```
with vector field ``v``, initial condition ``q_{0}`` and the solution
``q`` taking values in ``\\mathbb{R}^{d}``. Here, the vector field ``v``
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
    function v_i(t, q, v)
        v[1] = ...
        v[2] = ...
        ...
    end
```
where `t` is the current time, `q` is the current solution vector, and
`v` is the vector which holds the result of evaluating the vector field ``v_i``
on `t` and `q`.

"""
struct SODE{dType <: Number, tType <: Number, vType <: Tuple, N} <: Equation{dType, tType}
    d::Int
    n::Int
    v::vType
    t₀::tType
    q₀::Array{dType,N}
    periodicity::Vector{dType}

    function SODE{dType,tType,vType,N}(d, n, v, t₀, q₀; periodicity=[]) where {dType, tType, vType, N}
        @assert d == size(q₀,1)
        @assert n == size(q₀,2)
        @assert dType == eltype(q₀)
        @assert ndims(q₀) == N ∈ (1,2)

        if !(length(periodicity) == d)
            periodicity = zeros(dType, d)
        end

        new(d, n, v, t₀, q₀, periodicity)
    end
end

function SODE(v::VT, t₀::TT, q₀::DenseArray{DT}; periodicity=[]) where {VT,DT,TT}
    SODE{DT, TT, VT, ndims(q₀)}(size(q₀, 1), size(q₀, 2), v, t₀, q₀, periodicity=periodicity)
end

function SODE(v, q₀; periodicity=[])
    SODE(v, zero(eltype(q₀)), q₀, periodicity=periodicity)
end

Base.hash(SODE::SODE, h::UInt) = hash(SODE.d, hash(SODE.n, hash(SODE.v, hash(SODE.t₀, hash(SODE.q₀, hash(SODE.periodicity, h))))))
Base.:(==)(SODE1::SODE, SODE2::SODE) = (
                                SODE1.d == SODE2.d
                             && SODE1.n == SODE2.n
                             && SODE1.v == SODE2.v
                             && SODE1.t₀ == SODE2.t₀
                             && SODE1.q₀ == SODE2.q₀
                             && SODE1.periodicity == SODE2.periodicity)

function Base.similar(SODE::SODE{DT,TT,VT}, q₀::DenseArray{DT}) where {DT, TT, VT}
    similar(SODE, SODE.t₀, q₀)
end

function Base.similar(SODE::SODE{DT,TT,VT}, t₀::TT, q₀::DenseArray{DT}) where {DT, TT, VT}
    @assert SODE.d == size(q₀,1)
    SODE(SODE.v, t₀, q₀)
end
