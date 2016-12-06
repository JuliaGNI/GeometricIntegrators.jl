"""
`ODE`: Ordinary Differential Equation

Defines an initial value problem
```math
\\dot{q} (t) = f(t, q(t)) , \\qquad q(t_{0}) = q_{0} ,
```
with vector field ``f``, initial condition ``q_{0}`` and the solution
``q`` taking values in ``\\mathbb{R}^{d}``.

### Fields

* `d`: dimension of dynamical variable ``q`` and the vector field ``f``
* `f`: function computing the vector field
* `t₀`: initial time
* `q₀`: initial condition

The function `f` providing the vector field must have the interface
```julia
    function f(t, q, f)
        f[1] = ...
        f[2] = ...
        ...
    end
```
where `t` is the current time, `q` is the current solution vector, and
`f` is the vector which holds the result of evaluating the vector field ``f``
on `t` and `q`.

"""
immutable ODE{dType, tType, fType} <: Equation{dType, tType}
    d::Int
    n::Int
    f::fType
    t₀::tType
    q₀::Array{dType,2}

    function ODE(d, n, f, t₀, q₀)
        @assert d == size(q₀,1)
        @assert n == size(q₀,2)
        @assert dType == eltype(q₀)
        @assert ndims(q₀) ∈ (1,2)

        if ndims(q₀) == 1
            q₀ = reshape(q₀, d, n)
        end

        new(d, n, f, t₀, q₀)
    end
end

function ODE{DT,TT,FT}(f::FT, t₀::TT, q₀::DenseArray{DT})
    ODE{DT,TT,FT}(size(q₀, 1), size(q₀, 2), f, t₀, q₀)
end

function ODE(f, q₀)
    ODE(f, zero(eltype(q₀)), q₀)
end

Base.hash(ode::ODE, h::UInt) = hash(ode.d, hash(ode.n, hash(ode.f, hash(ode.t₀, hash(ode.q₀, h)))))
Base.:(==){DT1, DT2, TT1, TT2, FT1, FT2}(ode1::ODE{DT1,TT1,FT1}, ode2::ODE{DT2,TT2,FT2}) = (
                                ode1.d == ode2.d
                             && ode1.n == ode2.n
                             && ode1.f == ode2.f
                             && ode1.t₀ == ode2.t₀
                             && ode1.q₀ == ode2.q₀)
