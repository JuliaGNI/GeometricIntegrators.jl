
using ..CommonFunctions


struct LagrangePolynomial{T<:AbstractFloat}
    b::LagrangeBasis{T}
    c::Vector{T}

    function LagrangePolynomial{T}(b::LagrangeBasis, c) where {T}
        @assert b.n == length(c)
        new(b, c)
    end

    function LagrangePolynomial{T}(x::Vector, c) where {T}
        @assert length(x) == length(c)

        new(LagrangeBasis(x), c)
    end
end

function LagrangePolynomial(b::LagrangeBasis{T}, c::Vector{T}) where {T}
    LagrangePolynomial{T}(b, c)
end


function LagrangePolynomial(x::Vector{T}, c::Vector{T}) where {T}
    LagrangePolynomial{T}(x, c)
end


function Base.similar(lag::LagrangePolynomial{T}, c::Vector{T}) where {T}
    @assert length(c) == length(lag.c)
    LagrangePolynomial{T}(lag.b, c)
end


function CommonFunctions.evaluate!(pol::LagrangePolynomial{T}, x::Vector{T}, y::Vector{T}) where {T}
    @assert length(x) == length(y)

    local tx::T
    local ty::T

    for k in eachindex(x,y)
        tx = x[k]
        ty = 0
        for j in eachindex(pol.c)
            ty += pol.c[j] * evaluate(pol.b, j, tx)
        end
        y[k] = ty
    end

    return nothing
end

function CommonFunctions.evaluate!(b::LagrangeBasis{T}, c::Vector{T}, x::Vector{T}, y::Vector{T}) where {T}
    evaluate!(LagrangePolynomial(b, c), x, y)
end
