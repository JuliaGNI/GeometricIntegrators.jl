
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

function simd_scale!(x, a)
    @simd for i=1:length(x)
        @inbounds x[i] *= a
    end
end

"Copy the first dimension of a 2D array y into a 1D array x."
function simd_copy_xy_first!(x, y, j)
    @simd for i=1:length(x)
        @inbounds x[i] = y[i,j]
    end
end

"Copy a 1D array x into the first dimension of a 2D array y."
function simd_copy_yx_first!(x, y, j)
    @simd for i=1:length(x)
        @inbounds y[i,j] = x[i]
    end
end

function simd_axpy!(a, x, y)
    @simd for i=1:length(x)
        @inbounds y[i] += a*x[i]
    end
end

function simd_waxpy!(w, a, x, y)
    @simd for i=1:length(x)
        @inbounds w[i] = a*x[i] + y[i]
    end
end

function simd_mult!(w, X, y)
    @simd for i=1:length(w)
        w[i] = 0
        @simd for j=1:length(y)
            @inbounds w[i] += X[i,j] * y[j]
        end
    end
end
