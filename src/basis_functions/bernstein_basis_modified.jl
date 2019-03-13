
using OffsetArrays


struct BernsteinBasisModified{T,N} <: PolynomialBasis{T,N}
    x::OffsetArray{T,1,Array{T,1}}

    function BernsteinBasisModified{T,N}(x) where {T,N}
        @assert length(x) == N
        ox  = OffsetArray{T}(undef, 0:N-1)
        for i in eachindex(x)
            ox[i-1] = x[i]
        end
        new(ox)
    end
end

function BernsteinBasisModified(x::Vector{T}) where {T}
    BernsteinBasisModified{T,length(x)}(x)
end

nodes(b::BernsteinBasisModified{T,N})  where {T,N} = b.x.parent

Base.hash(b::BernsteinBasisModified, h::UInt) = hash(b.c, h)

Base.:(==)(b1::BernsteinBasisModified, b2::BernsteinBasisModified) = (b1.x == b2.x)

Base.isequal(b1::BernsteinBasisModified{T1,N1}, b2::BernsteinBasisModified{T2,N2}) where {T1,N1,T2,N2} = (b1 == b2 && T1 == T2 && N1 == N2)


function bernstein(b::BernsteinBasisModified{T,N}, i::Int, n::Int, x::T) where {T,N}
    if i < 0 || i > n
        return zero(T)
    elseif i == 0 && n == N-1
        (b.x[end] - x) / (b.x[end] - b.x[0])
    elseif i == n && n == N-1
        (x - b.x[0]) / (b.x[end] - b.x[0])
    else
        if n == 0
            return one(T)
        else
            return bernstein(b, i, n-1, x) * (1-x) + bernstein(b, i-1, n-1, x) * x
        end
    end
end


function eval_basis(b::BernsteinBasisModified{T,N}, i::Int, x::T) where {T,N}
    @assert i ≥ 1 && i ≤ N
    bernstein(b, i-1, N-1, x)
end


function deriv_basis(b::BernsteinBasisModified{T,N}, i::Int, x::T) where {T,N}
    if i == 1
        - 1 / (b.x[end] - b.x[0])
    elseif i == N
        + 1 / (b.x[end] - b.x[0])
    else
        @assert i ≥ 1 && i ≤ N
        (N-1) * ( bernstein(b, i-2, N-2, x) - bernstein(b, i-1, N-2, x) )
    end
end

deriv_basis(b::BernsteinBasisModified, i::Int, j::Int) = derivative(b, i, b.x[j])


# function int_basis(b::BernsteinBasisModified{T,N}, i::Int, x::T) where {T,N}
# end
#
# int_basis(b::BernsteinBasisModified, i::Int, j::Int) = integral(b, i, b.x[j])
