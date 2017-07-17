"""
`IODE`: Implicit Ordinary Differential Equation

Defines an implicit initial value problem
```math
\\begin{align*}
\\dot{q} (t) &= v(t) , &
q(t_{0}) &= q_{0} , \\\\
\\dot{p} (t) &= f(t, q(t), v(t)) , &
p(t_{0}) &= p_{0} , \\\\
p(t) &= α(t, q(t), v(t))
\\end{align*}
```
with vector field ``f``, the momentum defined by ``p``, initial conditions ``(q_{0}, p_{0})`` and the solution
``(q,p)`` taking values in ``\\mathbb{R}^{d} \\times \\mathbb{R}^{d}``.
This is a special case of a differential algebraic equation with dynamical
variables ``(q,p)`` and algebraic variable ``v``.

### Fields

* `d`: dimension of dynamical variables ``q`` and ``p`` as well as the vector fields ``f`` and ``p``
* `α`: function determining the momentum
* `f`: function computing the vector field
* `g`: function determining the projection, given by ∇α(q)λ
* `v`: function computing an initial guess for the velocity field (optional)
* `t₀`: initial time (optional)
* `q₀`: initial condition for `q`
* `p₀`: initial condition for `p`

The functions `α` and `f` must have the interface
```julia
    function α(t, q, v, p)
        p[1] = ...
        p[2] = ...
        ...
    end
```
and
```julia
    function f(t, q, v, f)
        f[1] = ...
        f[2] = ...
        ...
    end
```
where `t` is the current time, `q` is the current solution vector, `v` is the
current velocity and `f` and `p` are the vectors which hold the result of
evaluating the functions ``f`` and ``α`` on `t`, `q` and `v`.
The funtions `g` and `v` are specified by
```julia
    function g(t, q, λ, g)
        g[1] = ...
        g[2] = ...
        ...
    end
```
and
```julia
    function v(t, q, p, v)
        v[1] = ...
        v[2] = ...
        ...
    end
```
"""
struct IODE{dType <: Number, tType <: Number, αType <: Function, fType <: Function, gType <: Function, vType <: Function, N} <: Equation{dType, tType}
    d::Int
    m::Int
    n::Int
    α::αType
    f::fType
    g::gType
    v::vType
    t₀::tType
    q₀::Array{dType, N}
    p₀::Array{dType, N}
    periodicity::Vector{dType}

    function IODE{dType,tType,αType,fType,gType,vType,N}(d, n, α, f, g, v, t₀, q₀, p₀; periodicity=[]) where {dType <: Number, tType <: Number, αType <: Function, fType <: Function, gType <: Function, vType <: Function, N}
        @assert d == size(q₀,1) == size(p₀,1)
        @assert n == size(q₀,2) == size(p₀,2)
        @assert dType == eltype(q₀) == eltype(p₀)
        @assert ndims(q₀) == ndims(p₀) == N ∈ (1,2)

        if !(length(periodicity) == d)
            periodicity = zeros(dType, d)
        end

        new(d, d, n, α, f, g, v, t₀, q₀, p₀, periodicity)
    end
end


function IODE(α::ΑT, f::FT, g::GT, v::VT, t₀::TT, q₀::DenseArray{DT}, p₀::DenseArray{DT}; periodicity=[]) where {DT,TT,ΑT,FT,GT,VT}
    @assert size(q₀) == size(p₀)
    IODE{DT, TT, ΑT, FT, GT, VT, ndims(q₀)}(size(q₀, 1), size(q₀, 2), α, f, g, v, t₀, q₀, p₀, periodicity=periodicity)
end

function IODE(α::Function, f::Function, g::Function, t₀::Number, q₀::DenseArray, p₀::DenseArray; periodicity=[])
    IODE(α, f, g, function_v_dummy, t₀, q₀, p₀, periodicity=periodicity)
end

function IODE(α::Function, f::Function, g::Function, v::Function, q₀::DenseArray, p₀::DenseArray; periodicity=[])
    IODE(α, f, g, v, zero(eltype(q₀)), q₀, p₀, periodicity=periodicity)
end

function IODE(f, p, u, g, q₀, p₀; periodicity=[])
    IODE(α, f, g, function_v_dummy, zero(eltype(q₀)), q₀, p₀, periodicity=periodicity)
end

Base.hash(ode::IODE, h::UInt) = hash(ode.d, hash(ode.n, hash(ode.α, hash(ode.f, hash(ode.g, hash(ode.v, hash(ode.t₀, hash(ode.q₀, hash(ode.p₀, hash(ode.periodicity, h))))))))))
Base.:(==)(ode1::IODE, ode2::IODE) = (
                                ode1.d == ode2.d
                             && ode1.n == ode2.n
                             && ode1.α == ode2.α
                             && ode1.f == ode2.f
                             && ode1.g == ode2.g
                             && ode1.v == ode2.v
                             && ode1.t₀ == ode2.t₀
                             && ode1.q₀ == ode2.q₀
                             && ode1.p₀ == ode2.p₀
                             && ode1.periodicity == ode2.periodicity)
