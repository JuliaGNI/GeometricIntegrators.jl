"""
`IODE`: Implicit Ordinary Differential Equation

Defines an implicit initial value problem
```math
\\begin{align*}
\\dot{q} (t) &= v(t) , &
q(t_{0}) &= q_{0} , \\\\
\\dot{p} (t) &= f(t, q(t), v(t)) , &
p(t_{0}) &= p_{0} , \\\\
p(t) &= p(t, q(t), v(t))
\\end{align*}
```
with vector field ``f``, the momentum defined by ``p``, initial conditions ``(q_{0}, p_{0})`` and the solution
``(q,p)`` taking values in ``\\mathbb{R}^{d} \\times \\mathbb{R}^{d}``.
This is a special case of a differential algebraic equation with dynamical
variables ``(q,p)`` and algebraic variable ``v``.

### Fields

* `d`: dimension of dynamical variables ``q`` and ``p`` as well as the vector fields ``f`` and ``p``
* `f`: function computing the vector field
* `p`: function determining ``p``
* `t₀`: initial time
* `q₀`: initial condition for `q`
* `p₀`: initial condition for `p`

The functions `f` and `g` must have the interface
```julia
    function f(t, q, v, f)
        f[1] = ...
        f[2] = ...
        ...
    end
```
and
```julia
    function p(t, q, v, p)
        p[1] = ...
        p[2] = ...
        ...
    end
```
where `t` is the current time, `q` is the current solution vector, `v` is the
current velocity and `f` and `p` are the vectors which hold the result of
evaluating the functions ``f`` and ``p`` on `t`, `q` and `v`.
"""
immutable IODE{dType <: Number, tType <: Number, fType <: Function, pType <: Function, vType <: Function, N} <: Equation{dType, tType}
    d::Int
    n::Int
    f::fType
    p::pType
    v::vType
    t₀::tType
    q₀::Array{dType, N}
    p₀::Array{dType, N}

    function IODE(d, n, f, p, v, t₀, q₀, p₀)
        @assert d == size(q₀,1) == size(p₀,1)
        @assert n == size(q₀,2) == size(p₀,2)
        @assert dType == eltype(q₀) == eltype(p₀)
        @assert ndims(q₀) == ndims(p₀) == N ∈ (1,2)
        new(d, n, f, p, v, t₀, q₀, p₀)
    end
end


function IODE{DT,TT,FT,PT,VT}(f::FT, p::PT, v::VT, t₀::TT, q₀::DenseArray{DT}, p₀::DenseArray{DT})
    @assert size(q₀) == size(p₀)
    IODE{DT, TT, FT, PT, VT, ndims(q₀)}(size(q₀, 1), size(q₀, 2), f, p, v, t₀, q₀, p₀)
end

function IODE(f::Function, p::Function, t₀::Number, q₀::DenseArray, p₀::DenseArray)
    IODE(f, p, function_v_dummy, t₀, q₀, p₀)
end

function IODE(f::Function, p::Function, v::Function, q₀::DenseArray, p₀::DenseArray)
    IODE(f, p, v, zero(Float64), q₀, p₀)
end

function IODE(f, p, q₀, p₀)
    IODE(f, p, function_v_dummy, zero(Float64), q₀, p₀)
end

Base.hash(ode::IODE, h::UInt) = hash(ode.d, hash(ode.n, hash(ode.f, hash(ode.p, hash(ode.v, hash(ode.t₀, hash(ode.q₀, hash(ode.p₀, h))))))))
Base.:(==)(ode1::IODE, ode2::IODE) = (
                                ode1.d == ode2.d
                             && ode1.n == ode2.n
                             && ode1.f == ode2.f
                             && ode1.p == ode2.p
                             && ode1.v == ode2.v
                             && ode1.t₀ == ode2.t₀
                             && ode1.q₀ == ode2.q₀
                             && ode1.p₀ == ode2.p₀)
