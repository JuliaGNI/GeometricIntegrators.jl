
using Polynomials

using ..CommonFunctions


struct LegendreBasis{T, P, fT, pT, dT, iT} <: Basis{T,P}
    factors::fT
    polys::pT
    derivs::dT
    ints::iT
end

function LegendreBasis(T,P)

    factors = Tuple{}()
    polys   = Tuple{}()
    derivs  = Tuple{}()
    ints    = Tuple{}()

    for p in 0:P
        tfac, tpol = legendre_polynomial(p)
        tder = polyder(tpol)
        tint = polyint(tpol)

        factors = (factors..., tfac)
        polys   = (polys..., tpol)
        derivs  = (derivs..., tder)
        ints    = (ints..., tint)
    end


    LegendreBasis{T, P, typeof(factors), typeof(polys), typeof(derivs), typeof(ints)}(factors, polys, derivs, ints)
end

CommonFunctions.nbasis(b::LegendreBasis{T,P}) where {T,P} = P+1
CommonFunctions.nnodes(b::LegendreBasis{T,P}) where {T,P} = P+1
CommonFunctions.nodes(b::LegendreBasis{T,P})  where {T,P} = Vector{T}([])
CommonFunctions.degree(b::LegendreBasis{T,P}) where {T,P} = P

Base.hash(b::LegendreBasis{T,P}, h::UInt) where {T,P} = hash(T, hash(P, h))

Base.:(==)(b1::LegendreBasis{T1,P1}, b2::LegendreBasis{T2,P2}) where {T1,P1,T2,P2} = (T1 == T2 && P1 == P2)

Base.isequal(b1::LegendreBasis{T1,P1}, b2::LegendreBasis{T2,P2}) where {T1,P1,T2,P2} = (T1 == T2 && P1 == P2)


function legendre_polynomial(p, T=Float64)
    if p == 0
        (T(sqrt(1)), Poly([1]))
    elseif p == 1
        (T(sqrt(3)), Poly([-1, +2]))
    elseif p == 2
        (T(sqrt(5)), Poly([+1, -6, +6]))
    elseif p == 3
        (T(sqrt(7)), Poly([-1, +12, -30, +20]))
    elseif p == 4
        (T(sqrt(9)), Poly([+1, -20, +90, -140, +70]))
    elseif p == 5
        (T(sqrt(11)), Poly([-1, +30, -210, +560, -630, +252]))
    elseif p == 6
        (T(sqrt(13)), Poly([+1, -42, +420, -1680, +3150, -2772, +924]))
    elseif p == 7
        (T(sqrt(15)), Poly([-1, +56, -756, +4200, -11550, +16632, -12012, +3432]))
    elseif p == 8
        (T(sqrt(17)), Poly([+1, -72, +1260, -9240, +34650, -72072, +84084, -51480, +12870]))
    elseif p == 9
        (T(sqrt(19)), Poly([-1, +90, -1980, +18480, -90090, +252252, -420420, +411840, -218790, +48620]))
    elseif p == 10
        (T(sqrt(21)), Poly([+1, -110, +2970, -34320, +210210, -756756, +1681680, -2333760, +1969110, -923780, +184756]))
    else
        error("Legendre polynomials only implemented for 0 ≤ p ≤ 10.")
    end
end


function CommonFunctions.evaluate(b::LegendreBasis{T,P}, j::Int, x::T) where {T,P}
    return b.factors[j] * b.polys[j](x)
end


function derivative(b::LegendreBasis{T,P}, j::Int, x::T) where {T,P}
    return b.factors[j] * b.derivs[j](x)
end


function integral(b::LegendreBasis{T,P}, j::Int, x::T) where {T,P}
    return b.factors[j] * b.ints[j](x)
end
