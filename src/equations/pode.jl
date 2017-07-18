"""
`IODE`: Partitioned Ordinary Differential Equation

Defines a partitioned initial value problem
```math
\\begin{align*}
\\dot{q} (t) &= v(t, q(t), p(t)) , &
q(t_{0}) &= q_{0} , \\\\
\\dot{p} (t) &= f(t, q(t), p(t)) , &
p(t_{0}) &= p_{0} ,
\\end{align*}
```
with vector fields ``v`` and ``f``, initial conditions ``(q_{0}, p_{0})`` and the solution
``(q,p)`` taking values in ``\\mathbb{R}^{d} \\times \\mathbb{R}^{d}``.

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
struct PODE{dType <: Number, tType <: Number, vType <: Function, fType <: Function, N} <: Equation{dType, tType}
    d::Int
    n::Int
    v::vType
    f::fType
    t₀::tType
    q₀::Array{dType, N}
    p₀::Array{dType, N}
    periodicity::Vector{dType}

    function PODE{dType,tType,vType,fType,N}(d, n, v, f, t₀, q₀, p₀; periodicity=[]) where {dType <: Number, tType <: Number, vType <: Function, fType <: Function, N}
        @assert d == size(q₀,1) == size(p₀,1)
        @assert n == size(q₀,2) == size(p₀,2)
        @assert dType == eltype(q₀) == eltype(p₀)
        @assert ndims(q₀) == ndims(p₀) == N ∈ (1,2)

        if !(length(periodicity) == d)
            periodicity = zeros(dType, d)
        end

        new(d, n, v, f, t₀, q₀, p₀, periodicity)
    end
end


function PODE(v::VT, f::FT, t₀::TT, q₀::DenseArray{DT}, p₀::DenseArray{DT}; periodicity=[]) where {DT,TT,VT,FT}
    @assert size(q₀) == size(p₀)
    PODE{DT, TT, VT, FT, ndims(q₀)}(size(q₀, 1), size(q₀, 2), v, f, t₀, q₀, p₀, periodicity=periodicity)
end

function PODE(v, f, q₀, p₀; periodicity=[])
    PODE(v, f, zero(eltype(q₀)), q₀, p₀, periodicity=periodicity)
end

Base.hash(ode::PODE, h::UInt) = hash(ode.d, hash(ode.n, hash(ode.v, hash(ode.f, hash(ode.t₀, hash(ode.q₀, hash(ode.p₀, h)))))))
Base.:(==){DT1, DT2, TT1, TT2, VT1, VT2, FT1, FT2}(ode1::PODE{DT1,TT1,VT1,FT1}, ode2::PODE{DT2,TT2,VT2,FT2}) = (
                                ode1.d == ode2.d
                             && ode1.n == ode2.n
                             && ode1.v == ode2.v
                             && ode1.f == ode2.f
                             && ode1.t₀ == ode2.t₀
                             && ode1.q₀ == ode2.q₀
                             && ode1.p₀ == ode2.p₀
                             && ode1.periodicity == ode2.periodicity)
