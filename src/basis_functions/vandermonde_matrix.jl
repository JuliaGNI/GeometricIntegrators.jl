
function vandermonde_matrix{T}(x::Vector{T})
    local n = length(x)
	local V = Array{T,2}(n,n)

	for i in 1:n
		V[:,i] .= x.^(i-1)
	end

	V
end

function vandermonde_matrix_inverse{T}(x::Vector{T})
    local n = length(x)

    local L::Matrix{T} = zeros(n,n)
    local U::Matrix{T} = eye(n)
    local V::Matrix{T}

    L[1,1] = 1
    for i in 2:n
        for j in 1:i
            p = 1
            for k in 1:i
                if k ≠ j
                    p *= (x[j] - x[k])
                end
            end
            L[i,j] = 1/p
        end
    end

    i = 1
    for j in i+1:n
        U[i,j] = - U[i,j-1] * x[j-1]
    end

    for i in 2:n
        for j in i+1:n
            U[i,j] = U[i-1,j-1] - U[i,j-1] * x[j-1]
        end
    end

    V = *(U,L)
end
