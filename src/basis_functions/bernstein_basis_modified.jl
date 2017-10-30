
using OffsetArrays


struct BernsteinBasisModified{T,P} <: Basis{T,P}
    x::OffsetArray{T,1,Array{T,1}}

    function BernsteinBasisModified{T,P}(x) where {T,P}
        @assert length(x) == P+1
        ox  = OffsetArray(T, 0:P)
        for i in eachindex(x)
            ox[i-1] = x[i]
        end
        new(ox)
    end
end

function BernsteinBasisModified(x::Vector{T}) where {T}
    BernsteinBasisModified{T,length(x)-1}(x)
end

CommonFunctions.nbasis(b::BernsteinBasisModified{T,P}) where {T,P} = P+1
CommonFunctions.nnodes(b::BernsteinBasisModified{T,P}) where {T,P} = P+1
CommonFunctions.nodes(b::BernsteinBasisModified{T,P})  where {T,P} = b.x.parent
CommonFunctions.degree(b::BernsteinBasisModified{T,P}) where {T,P} = P

Base.hash(b::BernsteinBasisModified, h::UInt) = hash(b.c, h)

Base.:(==)(b1::BernsteinBasisModified, b2::BernsteinBasisModified) = (b1.x == b2.x)

Base.isequal(b1::BernsteinBasisModified{T1,P1}, b2::BernsteinBasisModified{T2,P2}) where {T1,P1,T2,P2} = (b1 == b2 && T1 == T2 && P1 == P2)


function bernstein(b::BernsteinBasisModified{T,P}, i::Int, n::Int, x::T) where {T,P}
    if i < 0 || i > n
        return zero(T)
    elseif i == 0 && n == P
        (b.x[end] - x) / (b.x[end] - b.x[0])
    elseif i == n && n == P
        (x - b.x[0]) / (b.x[end] - b.x[0])
    else
        if n == 0
            return one(T)
        else
            return bernstein(b, i, n-1, x) * (1-x) + bernstein(b, i-1, n-1, x) * x
        end
    end
end


function CommonFunctions.evaluate(b::BernsteinBasisModified{T,P}, i::Int, x::T) where {T,P}
    @assert i-1 ≥ 0 && i-1 ≤ P
    bernstein(b, i-1, P, x)
end


function CommonFunctions.derivative(b::BernsteinBasisModified{T,P}, i::Int, x::T) where {T,P}
    if i == 1
        - 1 / (b.x[end] - b.x[0])
    elseif i == P+1
        + 1 / (b.x[end] - b.x[0])
    else
        @assert i-1 ≥ 0 && i-1 ≤ P
        P * ( bernstein(b, i-2, P-1, x) - bernstein(b, i-1, P-1, x) )
    end
end

CommonFunctions.derivative(b::BernsteinBasisModified, i::Int, j::Int) = derivative(b, i, b.x[j])


# function CommonFunctions.integral(b::BernsteinBasisModified{T,P}, i::Int, x::T) where {T,P}
# end
#
# CommonFunctions.integral(b::BernsteinBasisModified, i::Int, j::Int) = integral(b, i, b.x[j])
