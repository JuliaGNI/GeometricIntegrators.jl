
using QuadratureRules

"Holds the coefficients of a projected Gauss-Legendre Runge-Kutta method."
struct CoefficientsPGLRK{T} <: AbstractCoefficients{T}
    @HeaderCoefficientsRK
    @CoefficientsRK

    P::Matrix{T}
    Q::Matrix{T}
    X::Matrix{T}
    W::Matrix{T}
    A::Matrix{T}

    function CoefficientsPGLRK{T}(name,o,s,a,b,c,P,X,W) where {T}
        @assert T <: Real
        @assert isa(name, Symbol)
        @assert isa(o, Integer)
        @assert isa(s, Integer)
        @assert s ≥ 2 "Number of stages must be ≥ 2"
        @assert s==size(a,1)==size(a,2)==length(b)==length(c)
        @assert s==size(P,1)==size(P,2)
        @assert s==size(X,1)==size(X,2)
        @assert s==size(W,1)==size(W,2)

        Q = inv(P)
        A = zero(a)
        B = zero(a)

        mul!(B, W, Q)
        mul!(A, P, B)

        new(name,o,s,a,b,c,P,Q,X,W,A)
    end
end

function CoefficientsPGLRK(name::Symbol, order::Int, a::Matrix{T}, b::Vector{T}, c::Vector{T}, P::Matrix{T}, X::Matrix{T}, W::Matrix{T}) where {T}
    CoefficientsPGLRK{T}(name, order, length(c), a, b, c, P, X, W)
end

function CoefficientsPGLRK(::Type{T}, s::Int) where {T}
    @assert s ≥ 2

    local ξᵢ::T

    # order
    o = 2s

    # obtain Gauss-Legendre nodes and weights
    gl_quad = GaussLegendreQuadrature(T,s)
    leg_basis = Legendre(T,s)

    a = zeros(T, s, s)
    t = zeros(T, s, s)
    P = zeros(T, s, s)
    X = zeros(T, s, s)
    W = zeros(T, s, s)

    for j in 1:s
        P[:,j] .= leg_basis[nodes(gl_quad), j-1]
    end

    Q = inv(P)

    X[1,1] = 0.5

    for i in 1:s-1
        ξᵢ = one(T) / sqrt(convert(T, 4i^2-1)) / 2
        X[i+1,i] = +ξᵢ
        X[i,i+1] = -ξᵢ
    end

    W[s, s-1] = +1
    W[s-1, s] = -1

    mul!(t, X, Q)
    mul!(a, P, t)

    CoefficientsPGLRK(Symbol("PGLRK", s), o, a, weights(gl_quad), nodes(gl_quad), P, X, W)
end

CoefficientsPGLRK(s::Int) = CoefficientsPGLRK(Float64, s)


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

# Print Runge-Kutta coefficients.
function Base.show(io::IO, tab::CoefficientsPGLRK)
    print(io, "Projected Gauss-Legendre Runge-Kutta Coefficients ", tab.name, " with ", tab.s, " stages and order ", tab.o)
    print(io, "  a = ", tab.a)
    print(io, "  b = ", tab.b)
    print(io, "  c = ", tab.c)
    print(io, "  P = ", tab.P)
    print(io, "  Q = ", tab.Q)
    print(io, "  X = ", tab.X)
    print(io, "  W = ", tab.W)
    print(io, "  A = ", tab.A)
end

function getTableauPGLRK(coeff::CoefficientsPGLRK{T}, λ, a::Matrix{T}) where {T}
    a .= coeff.a .+ λ .* coeff.A
end

function getTableauPGLRK(coeff::CoefficientsPGLRK{T}, λ) where {T}
    coeff.a .+ λ .* coeff.A
end
