
using ..CommonFunctions

immutable LagrangePolynomial{T<:AbstractFloat}
    p::Int
    n::Int
    x::Vector{T}
    y::Vector{T}
    denom::Vector{T}

    function LagrangePolynomial(x, y)
        local p::T

        @assert length(x) == length(y)

        denom = zeros(length(x))

        for i in 1:length(x)
            p = 1
            for j in 1:length(x)
                if i ≠ j
                    p *= (x[i] - x[j])
                end
            end
            denom[i] = 1/p
        end

        new(length(x)-1, length(x), x, y, denom)
    end

    function LagrangePolynomial(p, n, x, y, denom)
        @assert length(x) == length(y) == length(denom)
        new(p, n, x, y, denom)
    end
end

function LagrangePolynomial{T}(x::Vector{T}, y::Vector{T})
    LagrangePolynomial{T}(x, y)
end

function Base.similar{T}(lag::LagrangePolynomial{T}, y::Vector{T})
    @assert length(y) == lag.n
    LagrangePolynomial{T}(lag.p, lag.n, lag.x, y, lag.denom)
end

function CommonFunctions.evaluate!{T}(pol::LagrangePolynomial{T}, x::Vector{T}, y::Vector{T})
    @assert length(x) == length(y)

    local p::T
    local tx::T
    local ty::T

    for k in 1:length(x)
        tx = x[k]
        ty = 0
        for j in 1:pol.n
            p = 1
            for i in 1:pol.n
                if i ≠ j
                    p *= (tx - pol.x[i])
                end
            end
            ty += p * pol.y[j] * pol.denom[j]
        end
        y[k] = ty
    end

    return nothing
end

function CommonFunctions.evaluate!{T}(pol::LagrangePolynomial{T}, sy::Vector{T}, x::Vector{T}, y::Vector{T})
    spol = similar(pol, sy)
    evaluate!(spol, x, y)
end
