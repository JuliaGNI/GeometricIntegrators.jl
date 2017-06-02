
using ..CommonFunctions


struct LagrangeBasis{T<:AbstractFloat}
    p::Int
    n::Int
    x::Vector{T}

    denom::Vector{T}
    diffs::Matrix{T}

    function LagrangeBasis(x)
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

        new(n-1, n, x, denom)
    end

    function LagrangeBasis(p, n, x, denom)
        @assert length(x) == length(denom)
        new(p, n, x, denom)
    end
end

function LagrangeBasis{T}(x::Vector{T})
    LagrangeBasis{T}(x)
end


struct LagrangePolynomial{T<:AbstractFloat}
    b::LagrangeBasis{T}
    c::Vector{T}

    function LagrangePolynomial(b::LagrangeBasis, c)
        @assert b.n == length(c)
        new(b, c)
    end

    function LagrangePolynomial(x::Vector, c)
        @assert length(x) == length(c)

        new(LagrangeBasis(x), c)
    end
end

function LagrangePolynomial{T}(b::LagrangeBasis{T}, c::Vector{T})
    LagrangePolynomial{T}(b, c)
end


function LagrangePolynomial{T}(x::Vector{T}, c::Vector{T})
    LagrangePolynomial{T}(x, c)
end


function Base.similar{T}(lag::LagrangePolynomial{T}, c::Vector{T})
    @assert length(c) == length(lag.c)
    LagrangePolynomial{T}(lag.b, c)
end


function lagrange{T}(b::LagrangeBasis{T}, j::Int, x::T)
    local y::T = 1
    for i in 1:b.n
        i ≠ j ? y *= (x - b.x[i]) : nothing
    end
    y * b.denom[j]
end


function lagrange_derivative{T}(b::LagrangeBasis{T}, j::Int, x::T)
    local y::T = lagrange(b, j, x)
    for i in 1:b.n
        i ≠ j ? y /= (x - b.x[i]) : nothing
    end
    y
end


function lagrange_derivative{T}(b::LagrangeBasis{T}, j::Int)
    local y::T = 1
    for i in 1:b.n
        i ≠ j ? y /= b.diffs[j,i] : nothing
    end
    y
end


function CommonFunctions.evaluate!{T}(pol::LagrangePolynomial{T}, x::Vector{T}, y::Vector{T})
    @assert length(x) == length(y)

    local tx::T
    local ty::T

    for k in eachindex(x,y)
        tx = x[k]
        ty = 0
        for j in eachindex(pol.c)
            ty += pol.c[j] * lagrange(pol.b, j, tx)
        end
        y[k] = ty
    end

    return nothing
end

function CommonFunctions.evaluate!{T}(b::LagrangeBasis{T}, c::Vector{T}, x::Vector{T}, y::Vector{T})
    evaluate!(LagrangePolynomial(b, c), x, y)
end
