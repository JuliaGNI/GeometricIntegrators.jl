
using OffsetArrays

using ..CommonFunctions


struct BernsteinBasis{T,P} <: Basis{T,P}
    x::OffsetArray{T,1,Array{T,1}}

    function BernsteinBasis{T,P}(x) where {T,P}
        @assert length(x) == P+1
        ox  = OffsetArray(T, 0:P)
        for i in eachindex(x)
            ox[i-1] = x[i]
        end
        new(ox)
    end
end

function BernsteinBasis(x::Vector{T}) where {T}
    BernsteinBasis{T,length(x)-1}(x)
end

CommonFunctions.nbasis(b::BernsteinBasis{T,P}) where {T,P} = P+1
CommonFunctions.nnodes(b::BernsteinBasis{T,P}) where {T,P} = P+1
CommonFunctions.nodes(b::BernsteinBasis{T,P})  where {T,P} = b.x.parent
CommonFunctions.degree(b::BernsteinBasis{T,P}) where {T,P} = P

Base.hash(b::BernsteinBasis, h::UInt) = hash(b.c, h)

Base.:(==)(b1::BernsteinBasis, b2::BernsteinBasis) = (b1.x == b2.x)

Base.isequal(b1::BernsteinBasis{T1,P1}, b2::BernsteinBasis{T2,P2}) where {T1,P1,T2,P2} = (b1 == b2 && T1 == T2 && P1 == P2)


function bernstein(b::BernsteinBasis{T,P}, i::Int, n::Int, x::T) where {T,P}
    if i < 0 || i > n
        return zero(T)
    else
        if n == 0
            return one(T)
        else
            return bernstein(b, i, n-1, x) * (1-x) + bernstein(b, i-1, n-1, x) * x
        end
    end
end


function CommonFunctions.evaluate(b::BernsteinBasis{T,P}, i::Int, x::T) where {T,P}
    @assert i-1 ≥ 0 && i-1 ≤ P
    bernstein(b, i-1, P, x)
end


function derivative(b::BernsteinBasis{T,P}, i::Int, x::T) where {T,P}
    @assert i-1 ≥ 0 && i-1 ≤ P
    P * ( bernstein(b, i-2, P-1, x) - bernstein(b, i-1, P-1, x) )
end

derivative(b::BernsteinBasis, i::Int, j::Int) = derivative(b, i, b.x[j])


# function integral(b::BernsteinBasis{T,P}, i::Int, x::T) where {T,P}
# end
#
# integral(b::BernsteinBasis, i::Int, j::Int) = integral(b, i, b.x[j])
