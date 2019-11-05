
function compensated_summation(x::T, y::T, e::T) where {T}
    local err::T
    local res::T

    err = e   + x
    res = err + y
    e   = err + (y - res)

    return (res, e)
end

function compensated_summation(x::Vector{T}, y::Vector{T}, e::Vector{T}) where {T}
    local err::Vector{T}
    local res::Vector{T}
    local tmp::Vector{T}

    err = e   .+ x
    res = err .+ y
    tmp = y   .- res
    e  .= err .+ tmp

    return (res, e)
end


function L2norm(x)
    local l2::eltype(x) = 0
    for xᵢ in x
        l2 += xᵢ^2
    end
    l2
end

function l2norm(x)
    sqrt(L2norm(x))
end
