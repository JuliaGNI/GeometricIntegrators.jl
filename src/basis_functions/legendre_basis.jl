
import Polynomials 

using Polynomials: derivative, integrate


struct LegendreBasis{T, N, fT, pT, dT, iT} <: PolynomialBasis{T,N}
    factors::fT
    polys::pT
    derivs::dT
    ints::iT
end

function LegendreBasis(T,N)

    factors = Tuple{}()
    polys   = Tuple{}()
    derivs  = Tuple{}()
    ints    = Tuple{}()

    for p in 0:N-1
        tfac, tpol = legendre_polynomial(p)
        tder = derivative(tpol)
        tint = integrate(tpol)

        factors = (factors..., tfac)
        polys   = (polys..., tpol)
        derivs  = (derivs..., tder)
        ints    = (ints..., tint)
    end


    LegendreBasis{T, N, typeof(factors), typeof(polys), typeof(derivs), typeof(ints)}(factors, polys, derivs, ints)
end

nodes(b::LegendreBasis{T,N})  where {T,N} = Vector{T}([])

Base.hash(b::LegendreBasis{T,N}, h::UInt) where {T,N} = hash(T, hash(N, h))

Base.:(==)(b1::LegendreBasis{T1,N1}, b2::LegendreBasis{T2,N2}) where {T1,N1,T2,N2} = (T1 == T2 && N1 == N2)

Base.isequal(b1::LegendreBasis{T1,N1}, b2::LegendreBasis{T2,N2}) where {T1,N1,T2,N2} = (T1 == T2 && N1 == N2)


function legendre_polynomial(p, T=Float64)
    if p == 0
        (T(sqrt(1)), Polynomials.Polynomial([1]))
    elseif p == 1
        (T(sqrt(3)), Polynomials.Polynomial([-1, +2]))
    elseif p == 2
        (T(sqrt(5)), Polynomials.Polynomial([+1, -6, +6]))
    elseif p == 3
        (T(sqrt(7)), Polynomials.Polynomial([-1, +12, -30, +20]))
    elseif p == 4
        (T(sqrt(9)), Polynomials.Polynomial([+1, -20, +90, -140, +70]))
    elseif p == 5
        (T(sqrt(11)), Polynomials.Polynomial([-1, +30, -210, +560, -630, +252]))
    elseif p == 6
        (T(sqrt(13)), Polynomials.Polynomial([+1, -42, +420, -1680, +3150, -2772, +924]))
    elseif p == 7
        (T(sqrt(15)), Polynomials.Polynomial([-1, +56, -756, +4200, -11550, +16632, -12012, +3432]))
    elseif p == 8
        (T(sqrt(17)), Polynomials.Polynomial([+1, -72, +1260, -9240, +34650, -72072, +84084, -51480, +12870]))
    elseif p == 9
        (T(sqrt(19)), Polynomials.Polynomial([-1, +90, -1980, +18480, -90090, +252252, -420420, +411840, -218790, +48620]))
    elseif p == 10
        (T(sqrt(21)), Polynomials.Polynomial([+1, -110, +2970, -34320, +210210, -756756, +1681680, -2333760, +1969110, -923780, +184756]))
    else
        error("Legendre polynomials only implemented for 0 ≤ p ≤ 10.")
    end
end


function eval_basis(b::LegendreBasis{T,N}, j::Int, x::T) where {T,N}
    return b.factors[j] * b.polys[j](x)
end


function deriv_basis(b::LegendreBasis{T,N}, j::Int, x::T) where {T,N}
    return b.factors[j] * b.derivs[j](x)
end


function int_basis(b::LegendreBasis{T,N}, j::Int, x::T) where {T,N}
    return b.factors[j] * b.ints[j](x)
end
