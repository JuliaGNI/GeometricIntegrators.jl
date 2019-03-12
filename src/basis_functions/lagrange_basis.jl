
using Polynomials: Poly, polyint


struct LagrangeBasis{T,N} <: PolynomialBasis{T,N}
    x::Vector{T}

    denom::Vector{T}
    diffs::Matrix{T}
    vdminv::Matrix{T}

    function LagrangeBasis{T,N}(x) where {T,N}
        @assert length(x) == N

        local p::T

        denom = zeros(N)
        diffs = zeros(N,N)

        for i in 1:length(x)
            p = 1
            for j in 1:length(x)
                diffs[i,j] = x[i] - x[j]
                if i ≠ j
                    p *= diffs[i,j]
                end
            end
            denom[i] = 1/p
        end

        new(x, denom, diffs, vandermonde_matrix_inverse(x))
    end
end

function LagrangeBasis(x::Vector{T}) where {T}
    LagrangeBasis{T, length(x)}(x)
end

nodes(b::LagrangeBasis{T,N})  where {T,N} = b.x

Base.hash(b::LagrangeBasis, h::UInt) = hash(b.c, h)

Base.:(==)(b1::LagrangeBasis, b2::LagrangeBasis) = (b1.x == b2.x)

Base.isequal(b1::LagrangeBasis{T1,N1}, b2::LagrangeBasis{T2,N2}) where {T1,N1,T2,N2} = (T1 == T2 && N1 == N2)



function eval_basis(b::LagrangeBasis{T,N}, j::Int, x::T) where {T,N}
    local y::T = 1

    for i in 1:nnodes(b)
        i ≠ j ? y *= (x - b.x[i]) : nothing
    end
    y * b.denom[j]
end


function deriv_basis(b::LagrangeBasis{T,N}, j::Int, x::T) where {T,N}
    local y::T = 0
    local z::T

    for l in 1:nnodes(b)
        if l ≠ j
            z = 1 / b.diffs[j,l]
            for i in 1:nnodes(b)
                i ≠ j && i ≠ l ? z *= (x - b.x[i]) / b.diffs[j,i] : nothing
            end
            y += z
        end
    end
    y
end

deriv_basis(b::LagrangeBasis, j::Int, i::Int) = deriv_basis(b, j, b.x[i])


function int_basis(b::LagrangeBasis{T,N}, j::Int, x::T) where {T,N}
    y = zero(b.x)
    y[j] = 1
    lint = polyint(Poly(*(b.vdminv, y)))
    return lint(x)
end

int_basis(b::LagrangeBasis, j::Int, i::Int) = int_basis(b, j, b.x[i])
