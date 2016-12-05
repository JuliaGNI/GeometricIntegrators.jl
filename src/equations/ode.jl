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
    function f(t, q, fq)
        fq[1] = ...
        fq[2] = ...
        ...
    end
```
where `t` is the current time, `q` is the current solution vector, and
`fq` is the vector which holds the result of evaluating the vector field ``f``
on `t` and `q`.

"""
immutable ODE{T,F} <: Equation{T}
    d::Int
    n::Int
    f::F
    t₀::T
    q₀::Array{T,2}

    function ODE(d, n, f, t₀, q₀)
        @assert d == size(q₀,1)
        @assert n == size(q₀,2)
        @assert T == eltype(q₀)
        @assert ndims(q₀) ∈ (1,2)

        if ndims(q₀) == 1
            q₀ = reshape(q₀, d, n)
        end

        new(d, n, f, t₀, q₀)
    end
end

function ODE{T,F}(f::F, t₀::Real, q₀::DenseArray{T})
    ODE{T,F}(size(q₀, 1), size(q₀, 2), f, t₀, q₀)
end

function ODE(f, q₀)
    ODE(f, 0, q₀)
end

Base.hash(ode::ODE, h::UInt) = hash(ode.d, hash(ode.n, hash(ode.f, hash(ode.t₀, hash(ode.q₀, h)))))
Base.:(==){T1, T2, F1, F2}(ode1::ODE{T1,F1}, ode2::ODE{T2,F2}) = (
                                ode1.d == ode2.d
                             && ode1.n == ode2.n
                             && ode1.f == ode2.f
                             && ode1.t₀ == ode2.t₀
                             && ode1.q₀ == ode2.q₀)
