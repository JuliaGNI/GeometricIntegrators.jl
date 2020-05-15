
import Polynomials 

using Polynomials: derivative, integrate


struct LegendreBasisHierarchical{T, N, fT, pT, dT, iT} <: PolynomialBasis{T,N}
    x::Vector{T}
    factors::fT
    polys::pT
    derivs::dT
    ints::iT
end

function LegendreBasisHierarchical(T,N)

    if N < 2
        error("Number of Nodes of Hierarchical Legendre Basis must be at least 2 (Degree at least 1).")
    end

    x       = zeros(T,N)
    factors = Tuple{}()
    polys   = Tuple{}()
    derivs  = Tuple{}()
    ints    = Tuple{}()

    for p in 0:N-1
        tfac, tpol = legendre_polynomial_hierarchical(p, T)
        tder = derivative(tpol)
        tint = integrate(tpol)

        if p == 0
            x[p+1] = zero(T)
        elseif p == 1
            x[p+1] = one(T)
        else
            x[p+1] = T(0.5)
        end

        factors = (factors..., tfac)
        polys   = (polys..., tpol)
        derivs  = (derivs..., tder)
        ints    = (ints..., tint)
    end


    LegendreBasisHierarchical{T, N, typeof(factors), typeof(polys), typeof(derivs), typeof(ints)}(x, factors, polys, derivs, ints)
end

nodes(b::LegendreBasisHierarchical{T,N})  where {T,N} = b.x

Base.hash(b::LegendreBasisHierarchical{T,N}, h::UInt) where {T,N} = hash(T, hash(N, h))

Base.:(==)(b1::LegendreBasisHierarchical{T1,N1}, b2::LegendreBasisHierarchical{T2,N2}) where {T1,N1,T2,N2} = (T1 == T2 && N1 == N2)

Base.isequal(b1::LegendreBasisHierarchical{T1,N1}, b2::LegendreBasisHierarchical{T2,N2}) where {T1,N1,T2,N2} = (T1 == T2 && N1 == N2)


function legendre_polynomial_hierarchical(p, T=Float64)
    if p == 0
        (one(T), Polynomials.Polynomial([+1, -1]))
    elseif p == 1
        (one(T), Polynomials.Polynomial([ 0, +1]))
    else
        l1 = legendre_polynomial(p,   T)
        l2 = legendre_polynomial(p-2, T)

        l1p = l1[2]
        l2p = l2[2]

        (1/sqrt(2 * (2*p-1)), l1p - l2p)
    end
end


function eval_basis(b::LegendreBasisHierarchical{T,N}, j::Int, x::T) where {T,N}
    return b.factors[j] * b.polys[j](x)
end


function deriv_basis(b::LegendreBasisHierarchical{T,N}, j::Int, x::T) where {T,N}
    return b.factors[j] * b.derivs[j](x)
end


function int_basis(b::LegendreBasisHierarchical{T,N}, j::Int, x::T) where {T,N}
    return b.factors[j] * b.ints[j](x)
end
