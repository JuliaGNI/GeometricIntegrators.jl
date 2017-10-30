
using Polynomials


struct LagrangeBasis{T,P} <: Basis{T,P}
    x::Vector{T}

    denom::Vector{T}
    diffs::Matrix{T}
    vdminv::Matrix{T}

    function LagrangeBasis{T,P}(x) where {T,P}
        @assert length(x) == P+1

        local p::T

        denom = zeros(P+1)
        diffs = zeros(P+1, P+1)

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
    LagrangeBasis{T, length(x)-1}(x)
end

CommonFunctions.nbasis(b::LagrangeBasis{T,P}) where {T,P} = P+1
CommonFunctions.nnodes(b::LagrangeBasis{T,P}) where {T,P} = P+1
CommonFunctions.nodes(b::LagrangeBasis{T,P})  where {T,P} = b.x
CommonFunctions.degree(b::LagrangeBasis{T,P}) where {T,P} = P

Base.hash(b::LagrangeBasis, h::UInt) = hash(b.c, h)

Base.:(==)(b1::LagrangeBasis, b2::LagrangeBasis) = (b1.x == b2.x)

Base.isequal(b1::LagrangeBasis{T1,P1}, b2::LagrangeBasis{T2,P2}) where {T1,P1,T2,P2} = (T1 == T2 && P1 == P2)


function CommonFunctions.evaluate(b::LagrangeBasis{T,P}, j::Int, x::T) where {T,P}
    local y::T = 1

    for i in 1:nnodes(b)
        i ≠ j ? y *= (x - b.x[i]) : nothing
    end
    y * b.denom[j]
end


function CommonFunctions.derivative(b::LagrangeBasis{T,P}, j::Int, x::T) where {T,P}
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

CommonFunctions.derivative(b::LagrangeBasis, j::Int, i::Int) = derivative(b, j, b.x[i])


function CommonFunctions.integral(b::LagrangeBasis{T,P}, j::Int, x::T) where {T,P}
    y = zeros(b.x)
    y[j] = 1
    lint = polyint(Poly(*(b.vdminv, y)))
    return lint(x)
end

CommonFunctions.integral(b::LagrangeBasis, j::Int, i::Int) = integral(b, j, b.x[i])
