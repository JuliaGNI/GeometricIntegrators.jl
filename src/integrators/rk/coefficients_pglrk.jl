@doc raw"""
Holds the coefficients of a projected Gauss-Legendre Runge-Kutta method.

The construction follows the W-transformation: with the Legendre basis
``\hat{P}_{j}`` orthonormal on ``[0,1]`` evaluated at the ``s`` Gauss-Legendre nodes
``c_i``,
```math
P_{ij} = \hat{P}_{j-1} (c_i) ,
\qquad
Q = P^{-1} ,
```
the Gauss tableau itself is recovered as ``a = P X Q`` from the tridiagonal
```math
X_{11} = \tfrac{1}{2} ,
\qquad
X_{i+1,i} = + \xi_i ,
\qquad
X_{i,i+1} = - \xi_i ,
\qquad
\xi_i = \frac{1}{2 \sqrt{4 i^2 - 1}} ,
```
and a skew perturbation ``A = P W Q`` is built from the skew ``W`` with the single
nonzero pair ``W_{s,s-1} = +1``, ``W_{s-1,s} = -1``. This yields the one-parameter
family of tableaus
```math
a(\lambda) = a + \lambda A ,
```
used by `PGLRK` (with ``\lambda`` chosen to conserve energy) and by `VPRKpTableau`
(with ``\lambda`` chosen to enforce the Dirac constraint).

Two structural properties of ``A`` matter, and both are asserted in the tests:

* ``B A`` is **skew** for ``B = \mathrm{diag}(b)``, because ``P^T B P = I``. Hence
  ``b_i A_{ij} + b_j A_{ji} = 0`` and ``a(\lambda)`` satisfies the symplecticity
  condition for **every** ``\lambda`` and every choice of nonzero slot in ``W``.
* ``b^T A = 0`` and ``A \mathbb{1} = 0``, so ``a(\lambda)`` retains the quadrature and
  row-sum conditions — and therefore the order — but this holds **only** for
  ``s \geq 3``. At ``s = 2`` the pair ``W_{2,1} / W_{1,2}`` coincides with the
  order-determining ``\xi_1`` entries of ``X``, and the perturbation destroys
  consistency. The ``(2,1)`` slot is the only offending slot, at every ``s``.

For that reason ``s \geq 3`` is required. A two-stage method would remain symplectic
but lose consistency, i.e. it would be a symplectic method of reduced order.

### Reference

Luigi Brugnano, Felice Iavernaro and Donato Trigiante.
Energy- and quadratic invariants-preserving integrators based upon Gauss collocation
formulae.
SIAM Journal on Numerical Analysis 50(6), pp. 2897-2916, 2012.
doi: [10.1137/110856617](https://doi.org/10.1137/110856617).
"""
struct CoefficientsPGLRK{T} <: AbstractCoefficients{T}
    @HeaderCoefficientsRK
    @CoefficientsRK

    P::Matrix{T}
    Q::Matrix{T}
    X::Matrix{T}
    W::Matrix{T}
    A::Matrix{T}

    function CoefficientsPGLRK{T}(name, o, s, a, b, c, P, X, W) where {T}
        @assert T <: Real
        @assert isa(name, Symbol)
        @assert isa(o, Integer)
        @assert isa(s, Integer)
        @assert s ≥ 3 "Number of stages must be ≥ 3: the skew perturbation of a " *
                      "two-stage Gauss tableau occupies the order-determining " *
                      "(2,1)/(1,2) entries of X and destroys consistency."
        @assert s == size(a, 1) == size(a, 2) == length(b) == length(c)
        @assert s == size(P, 1) == size(P, 2)
        @assert s == size(X, 1) == size(X, 2)
        @assert s == size(W, 1) == size(W, 2)

        Q = inv(P)
        A = zero(a)
        B = zero(a)

        mul!(B, W, Q)
        mul!(A, P, B)

        new(name, o, s, a, b, c, zero(a), zero(b), zero(c), P, Q, X, W, A)
    end
end

function CoefficientsPGLRK(name::Symbol, order::Int, a::Matrix{T}, b::Vector{T}, c::Vector{T}, P::Matrix{T}, X::Matrix{T}, W::Matrix{T}) where {T}
    CoefficientsPGLRK{T}(name, order, length(c), a, b, c, P, X, W)
end

function CoefficientsPGLRK(::Type{T}, s::Int) where {T}
    # Checked here as well as in the inner constructor, because `W[s,s-1]` below
    # would raise a BoundsError for s < 2 before the inner assertion is reached.
    @assert s ≥ 3 "Number of stages must be ≥ 3: the skew perturbation of a " *
                  "two-stage Gauss tableau occupies the order-determining " *
                  "(2,1)/(1,2) entries of X and destroys consistency."

    local ξᵢ::T

    # order
    o = 2s

    # obtain Gauss-Legendre nodes and weights and the Legendre basis
    gl_quad = GaussLegendreQuadrature(T, s)
    gl_nodes = QuadratureRules.nodes(gl_quad)
    gl_weights = QuadratureRules.weights(gl_quad)
    leg_basis = Legendre(T, s)

    a = zeros(T, s, s)
    t = zeros(T, s, s)
    P = zeros(T, s, s)
    X = zeros(T, s, s)
    W = zeros(T, s, s)

    # Legendre-Vandermonde matrix of the W-transformation.
    # `Legendre` is indexed from zero, hence `j-1`.
    for j in 1:s
        P[:, j] .= leg_basis[gl_nodes, j-1]
    end

    X[1, 1] = one(T) / 2

    for i in 1:s-1
        ξᵢ = one(T) / sqrt(convert(T, 4i^2 - 1)) / 2
        X[i+1, i] = +ξᵢ
        X[i, i+1] = -ξᵢ
    end

    W[s, s-1] = +1
    W[s-1, s] = -1

    mul!(t, X, inv(P))
    mul!(a, P, t)

    CoefficientsPGLRK(Symbol("PGLRK", s), o, a, gl_weights, gl_nodes, P, X, W)
end

CoefficientsPGLRK(s::Int) = CoefficientsPGLRK(Float64, s)


GeometricBase.order(coeff::CoefficientsPGLRK) = coeff.o
GeometricBase.description(::CoefficientsPGLRK) = "Projected Gauss-Legendre Runge-Kutta coefficients"

RungeKutta.nstages(coeff::CoefficientsPGLRK) = coeff.s
RungeKutta.eachstage(coeff::CoefficientsPGLRK) = 1:coeff.s


Base.hash(tab::CoefficientsPGLRK, h::UInt) = hash(tab.o, hash(tab.s, hash(tab.a, hash(tab.b, hash(tab.c, hash(tab.P, hash(tab.Q, hash(tab.X, hash(tab.W, hash(tab.A, h))))))))))

Base.:(==)(tab1::CoefficientsPGLRK, tab2::CoefficientsPGLRK) = (tab1.o == tab2.o
                                                             && tab1.s == tab2.s
                                                             && tab1.a == tab2.a
                                                             && tab1.b == tab2.b
                                                             && tab1.c == tab2.c
                                                             && tab1.P == tab2.P
                                                             && tab1.Q == tab2.Q
                                                             && tab1.X == tab2.X
                                                             && tab1.W == tab2.W
                                                             && tab1.A == tab2.A)

Base.isequal(tab1::CoefficientsPGLRK{T1}, tab2::CoefficientsPGLRK{T2}) where {T1,T2} = (tab1 == tab2 && T1 == T2 && typeof(tab1) == typeof(tab2))

Base.isapprox(tab1::CoefficientsPGLRK, tab2::CoefficientsPGLRK; kwargs...) = (tab1.o == tab2.o
                                                             && tab1.s == tab2.s
                                                             && isapprox(tab1.a, tab2.a; kwargs...)
                                                             && isapprox(tab1.b, tab2.b; kwargs...)
                                                             && isapprox(tab1.c, tab2.c; kwargs...)
                                                             && isapprox(tab1.A, tab2.A; kwargs...))

# Print Runge-Kutta coefficients.
function Base.show(io::IO, tab::CoefficientsPGLRK)
    print(io, "Projected Gauss-Legendre Runge-Kutta Coefficients ", tab.name, " with ", tab.s, " stages and order ", tab.o, "\n")
    print(io, "  a = ", tab.a, "\n")
    print(io, "  b = ", tab.b, "\n")
    print(io, "  c = ", tab.c, "\n")
    print(io, "  P = ", tab.P, "\n")
    print(io, "  Q = ", tab.Q, "\n")
    print(io, "  X = ", tab.X, "\n")
    print(io, "  W = ", tab.W, "\n")
    print(io, "  A = ", tab.A, "\n")
end


@doc raw"""
```julia
getTableauPGLRK(coeff::CoefficientsPGLRK, λ, a::AbstractMatrix)
getTableauPGLRK(coeff::CoefficientsPGLRK, λ)
```

Evaluate the one-parameter family of tableaus ``a(\lambda) = a + \lambda A`` of a
[`CoefficientsPGLRK`](@ref). The three-argument form writes into `a` in place and
accepts any `AbstractMatrix` destination, so that `λ` may be a dual number during
automatic differentiation.
"""
function getTableauPGLRK(coeff::CoefficientsPGLRK, λ, a::AbstractMatrix)
    a .= coeff.a .+ λ .* coeff.A
end

function getTableauPGLRK(coeff::CoefficientsPGLRK, λ)
    coeff.a .+ λ .* coeff.A
end
