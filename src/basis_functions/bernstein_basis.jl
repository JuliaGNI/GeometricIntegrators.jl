
using OffsetArrays


struct BernsteinBasis{T,N} <: PolynomialBasis{T,N}
    x::OffsetArray{T,1,Array{T,1}}

    function BernsteinBasis{T,N}(x) where {T,N}
        @assert length(x) == N
        ox  = OffsetArray(T, 0:N-1)
        for i in eachindex(x)
            ox[i-1] = x[i]
        end
        new(ox)
    end
end

function BernsteinBasis(x::Vector{T}) where {T}
    BernsteinBasis{T,length(x)}(x)
end

nodes(b::BernsteinBasis{T,N})  where {T,N} = b.x.parent

Base.hash(b::BernsteinBasis, h::UInt) = hash(b.c, h)

Base.:(==)(b1::BernsteinBasis, b2::BernsteinBasis) = (b1.x == b2.x)

Base.isequal(b1::BernsteinBasis{T1,N1}, b2::BernsteinBasis{T2,N2}) where {T1,N1,T2,N2} = (b1 == b2 && T1 == T2 && N1 == N2)


function bernstein(b::BernsteinBasis{T,N}, i::Int, n::Int, x::T) where {T,N}
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


function eval_basis(b::BernsteinBasis{T,N}, i::Int, x::T) where {T,N}
    @assert i ≥ 1 && i ≤ N
    bernstein(b, i-1, N, x)
end


function deriv_basis(b::BernsteinBasis{T,N}, i::Int, x::T) where {T,N}
    @assert i ≥ 1 && i ≤ N
    P * ( bernstein(b, i-2, N-1, x) - bernstein(b, i-1, N-1, x) )
end

deriv_basis(b::BernsteinBasis, i::Int, j::Int) = derivative(b, i, b.x[j])


# function int_basis(b::BernsteinBasis{T,N}, i::Int, x::T) where {T,N}
# end
#
# int_basis(b::BernsteinBasis, i::Int, j::Int) = integral(b, i, b.x[j])
