
function istriustrict(A::Matrix{T}) where {T}
    m, n = size(A)
    @inbounds for j = 1:min(n,m-1), i = j:m
        if A[i,j] ≠ zero(T)
            return false
        end
    end
    return true
end

function istrilstrict(A::Matrix{T}) where {T}
    m, n = size(A)
    @inbounds for j = 2:n, i = 1:min(j,m)
        if A[i,j] ≠ zero(T)
            return false
        end
    end
    return true
end


function symplectic_matrix(t, q::AbstractVector{DT}, p::AbstractVector{DT}) where {DT}
    @assert length(q) == length(p)

    local d = length(q)

    ω = zeros(DT, 2d, 2d)

    for i in 1:d
        ω[i, d+i] = +1
        ω[d+i, i] = -1
    end

    return ω
end
