
using Polynomials

using ..CommonFunctions


struct LagrangeBasis{T<:AbstractFloat}
    p::Int
    n::Int
    x::Vector{T}

    denom::Vector{T}
    diffs::Matrix{T}
    vdminv::Matrix{T}

    function LagrangeBasis{T}(x) where {T}
        local p::T
        local n = length(x)

        denom = zeros(n)
        diffs = zeros(n,n)

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

        new(n-1, n, x, denom, diffs, vandermonde_matrix_inverse(x))
    end
end

function LagrangeBasis(x::Vector{T}) where {T}
    LagrangeBasis{T}(x)
end


Base.hash(b::LagrangeBasis, h::UInt) = hash(b.p, hash(b.n, hash(b.c, h)))

Base.:(==)(b1::LagrangeBasis, b2::LagrangeBasis) = (b1.p == b2.p
                                                 && b1.n == b2.n
                                                 && b1.x == b2.x)

Base.isequal(b1::LagrangeBasis{T1}, b2::LagrangeBasis{T2}) where {T1,T2} = (b1 == b2 && T1 == T2)


function CommonFunctions.evaluate(b::LagrangeBasis{T}, j::Int, x::T) where {T}
    local y::T = 1

    for i in 1:b.n
        i ≠ j ? y *= (x - b.x[i]) : nothing
    end
    y * b.denom[j]
end


function derivative(b::LagrangeBasis{T}, j::Int, x::T) where {T}
    local y::T = 0
    local z::T

    for l in 1:b.n
        if l ≠ j
            z = b.diffs[j,l]
            for i in 1:b.n
                i ≠ j && i ≠ l ? z *= (x - b.x[i]) / b.diffs[j,i] : nothing
            end
            y += z
        end
    end
    y
end

derivative(b::LagrangeBasis, j::Int, i::Int) = derivative(b, j, b.x[i])


function integral(b::LagrangeBasis{T}, j::Int, x::T) where {T}
    y = zeros(b.x)
    y[j] = 1
    lint = polyint(Poly(*(b.vdminv, y)))
    return lint(x)
end

integral(b::LagrangeBasis, j::Int, i::Int) = integral(b, j, b.x[i])
