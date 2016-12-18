"""
`IODE`: Implicit Ordinary Differential Equation

Defines an implicit initial value problem
```math
\\begin{align*}
\\dot{q} (t) &= v(t) , &
q(t_{0}) &= q_{0} , \\\\
\\dot{p} (t) &= f(t, q(t), v(t)) , &
p(t_{0}) &= p_{0} , \\\\
p(t) &= g(t, q(t), v(t))
\\end{align*}
```
with vector fields ``f`` and ``g``, initial conditions ``(q_{0}, p_{0})`` and the solution
``(q,p)`` taking values in ``\\mathbb{R}^{d} \\times \\mathbb{R}^{d}``.
This is a special case of a differential algebraic equation with dynamical
variables ``(q,p)`` and algebraic variable ``v``.

### Fields

* `d`: dimension of dynamical variables ``q`` and ``p`` as well as the vector fields ``f`` and ``g``
* `f`: function computing the vector field
* `g`: function determining the algebraic constraint
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
    function g(t, q, v, g)
        g[1] = ...
        g[2] = ...
        ...
    end
```
where `t` is the current time, `q` is the current solution vector, `v` is the
current velocity and `f` and `g` are the vectors which hold the result of
evaluating the functions ``f`` and ``g`` on `t`, `q` and `v`.
"""
immutable IODE{dType <: Number, tType <: Number, fType <: Function, gType <: Function} <: Equation{dType, tType}
    d::Int
    n::Int
    f::fType
    g::gType
    t₀::tType
    q₀::Array{dType, 2}
    p₀::Array{dType, 2}

    function IODE(d, n, f, g, t₀, q₀, p₀)
        @assert d == size(q₀,1) == size(p₀,1)
        @assert n == size(q₀,2) == size(p₀,2)

        @assert dType == eltype(q₀) == eltype(p₀)

        @assert ndims(q₀) ∈ (1,2)
        @assert ndims(p₀) ∈ (1,2)

        if ndims(q₀) == 1
            q₀ = reshape(q₀, d, n)
        end

        if ndims(p₀) == 1
            p₀ = reshape(p₀, d, n)
        end

        new(d, n, f, g, t₀, q₀, p₀)
    end
end


function IODE{DT,TT,FT,GT}(f::FT, g::GT, t₀::TT, q₀::DenseArray{DT}, p₀::DenseArray{DT})
    @assert size(q₀) == size(p₀)
    IODE{DT,TT,FT,GT}(size(q₀, 1), size(q₀, 2), f, g, t₀, q₀, p₀)
end

function IODE(f, g, q₀, p₀)
    IODE(f, g, zero(eltype(q₀)), q₀, p₀)
end

Base.hash(ode::IODE, h::UInt) = hash(ode.d, hash(ode.n, hash(ode.f, hash(ode.g, hash(ode.t₀, hash(ode.q₀, hash(ode.p₀, h)))))))
Base.:(==)(ode1::IODE, ode2::IODE) = (
                                ode1.d == ode2.d
                             && ode1.n == ode2.n
                             && ode1.f == ode2.f
                             && ode1.g == ode2.g
                             && ode1.t₀ == ode2.t₀
                             && ode1.q₀ == ode2.q₀
                             && ode1.p₀ == ode2.p₀)
