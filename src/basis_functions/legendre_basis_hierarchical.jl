
using Polynomials

using ..CommonFunctions


struct LegendreBasisHierarchical{T, P, fT, pT, dT, iT} <: Basis{T,P}
    x::Vector{T}
    factors::fT
    polys::pT
    derivs::dT
    ints::iT
end

function LegendreBasisHierarchical(T,P)

    if P < 1
        error("Degree of Hierarchical Legendre Basis must be at least 1.")
    end

    x       = zeros(T,P+1)
    factors = Tuple{}()
    polys   = Tuple{}()
    derivs  = Tuple{}()
    ints    = Tuple{}()

    for p in 0:P
        tfac, tpol = legendre_polynomial_hierarchical(p, T)
        tder = polyder(tpol)
        tint = polyint(tpol)

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


    LegendreBasisHierarchical{T, P, typeof(factors), typeof(polys), typeof(derivs), typeof(ints)}(x, factors, polys, derivs, ints)
end

CommonFunctions.nbasis(b::LegendreBasisHierarchical{T,P}) where {T,P} = P+1
CommonFunctions.nnodes(b::LegendreBasisHierarchical{T,P}) where {T,P} = P+1
CommonFunctions.nodes(b::LegendreBasisHierarchical{T,P})  where {T,P} = b.x
CommonFunctions.degree(b::LegendreBasisHierarchical{T,P}) where {T,P} = P

Base.hash(b::LegendreBasisHierarchical{T,P}, h::UInt) where {T,P} = hash(T, hash(P, h))

Base.:(==)(b1::LegendreBasisHierarchical{T1,P1}, b2::LegendreBasisHierarchical{T2,P2}) where {T1,P1,T2,P2} = (T1 == T2 && P1 == P2)

Base.isequal(b1::LegendreBasisHierarchical{T1,P1}, b2::LegendreBasisHierarchical{T2,P2}) where {T1,P1,T2,P2} = (T1 == T2 && P1 == P2)


function legendre_polynomial_hierarchical(p, T=Float64)
    if p == 0
        (one(T), Poly([+1, -1]))
    elseif p == 1
        (one(T), Poly([ 0, +1]))
    else
        l1 = legendre_polynomial(p,   T)
        l2 = legendre_polynomial(p-2, T)

        l1p = l1[2]
        l2p = l2[2]

        (1/sqrt(2 * (2*p-1)), l1p - l2p)
    end
end


function CommonFunctions.evaluate(b::LegendreBasisHierarchical{T,P}, j::Int, x::T) where {T,P}
    return b.factors[j] * b.polys[j](x)
end


function derivative(b::LegendreBasisHierarchical{T,P}, j::Int, x::T) where {T,P}
    return b.factors[j] * b.derivs[j](x)
end


function integral(b::LegendreBasisHierarchical{T,P}, j::Int, x::T) where {T,P}
    return b.factors[j] * b.ints[j](x)
end
