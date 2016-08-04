
function istriustrict(A::AbstractMatrix)
    m, n = size(A)
    for j = 1:min(n,m-1), i = j:m
        if A[i,j] != 0
            return false
        end
    end
    return true
end

function istrilstrict(A::AbstractMatrix)
    m, n = size(A)
    for j = 2:n, i = 1:min(j,m)
        if A[i,j] != 0
            return false
        end
    end
    return true
end
