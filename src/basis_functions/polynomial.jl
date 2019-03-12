
struct Polynomial{DT<:AbstractFloat, BT <: Basis}
    b::BT
    c::Vector{DT}

    function Polynomial{DT,BT}(b::BT, c) where {DT, BT <: Basis{DT}}
        @assert nnodes(b) == length(c)
        new(b, c)
    end
end

function Polynomial(b::Basis{DT}, c::Vector{DT}) where {DT}
    Polynomial{DT, typeof(b)}(b, c)
end


function LagrangePolynomial(x::Vector{DT}, c::Vector{DT}) where {DT}
    b = LagrangeBasis(x)
    Polynomial{DT, typeof(b)}(b, c)
end


function BernsteinPolynomial(x::Vector{DT}, c::Vector{DT}) where {DT}
    b = BernsteinBasis(x)
    Polynomial{DT, typeof(b)}(b, c)
end


function Base.similar(lag::Polynomial{T}, c::Vector{T}) where {T}
    @assert length(c) == length(lag.c)
    Polynomial(lag.b, c)
end


function evaluate!(pol::Polynomial{T}, x::Vector{T}, y::Vector{T}) where {T}
    @assert length(x) == length(y)

    local tx::T
    local ty::T

    for k in eachindex(x,y)
        tx = x[k]
        ty = 0
        for j in eachindex(pol.c)
            ty += pol.c[j] * eval_basis(pol.b, j, tx)
        end
        y[k] = ty
    end

    return nothing
end

function evaluate!(b::LagrangeBasis{T}, c::Vector{T}, x::Vector{T}, y::Vector{T}) where {T}
    evaluate!(Polynomial(b, c), x, y)
end
