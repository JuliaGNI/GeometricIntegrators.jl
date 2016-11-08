
function istriustrict(A::Matrix)
    m, n = size(A)
    for j = 1:min(n,m-1), i = j:m
        if A[i,j] != 0.
            return false
        end
    end
    return true
end

function istrilstrict(A::Matrix)
    m, n = size(A)
    for j = 2:n, i = 1:min(j,m)
        if A[i,j] != 0.
            return false
        end
    end
    return true
end

function simd_scale!(x, a)
    @simd for i=1:length(x)
        @inbounds x[i] *= a
    end
    nothing
end

function simd_copy!(x, y)
    @assert length(x) == length(y)
    @simd for i=1:length(y)
        @inbounds y[i] = x[i]
    end
    nothing
end

"Copy the first dimension of a 2D array y into a 1D array x."
function simd_copy_xy_first!(x, y, j)
    @assert length(x) == size(y, 1)
    @simd for i=1:length(x)
        @inbounds x[i] = y[i,j]
    end
    nothing
end

"Copy the first dimension of a 3D array y into a 1D array x."
function simd_copy_xy_first!(x, y, j, k)
    @assert length(x) == size(y, 1)
    @simd for i=1:length(x)
        @inbounds x[i] = y[i,j,k]
    end
    nothing
end

"Copy a 1D array x into the first dimension of a 2D array y."
function simd_copy_yx_first!(x, y, j)
    @assert length(x) == size(y, 1)
    @simd for i=1:size(y, 1)
        @inbounds y[i,j] = x[i]
    end
    nothing
end

"Copy a 1D array x into the first dimension of a 3D array y."
function simd_copy_yx_first!(x, y, j, k)
    @assert length(x) == size(y, 1)
    @simd for i=1:size(y, 1)
        @inbounds y[i,j,k] = x[i]
    end
    nothing
end

function simd_xpy!(x, y)
    @assert length(x) == length(y)
    @simd for i=1:length(y)
        @inbounds y[i] += x[i]
    end
    nothing
end

function simd_axpy!(a, x, y)
    @assert length(x) == length(y)
    @simd for i=1:length(y)
        @inbounds y[i] += a*x[i]
    end
    nothing
end

function simd_wxpy!(w, x, y)
    @assert length(x) == length(y) == length(w)
    @simd for i=1:length(w)
        @inbounds w[i] = x[i] + y[i]
    end
    nothing
end

function simd_waxpy!(w, a, x, y)
    @assert length(x) == length(y) == length(w)
    simd_copy!(y, w)
    simd_axpy!(a, x, w)
    nothing
end

function simd_mult!(w, X, y)
    @assert length(w) == size(X, 1)
    @assert length(y) == size(X, 2)
    for i=1:length(w)
        w[i] = 0.
        @simd for j=1:length(y)
            @inbounds w[i] += X[i,j] * y[j]
        end
    end
    nothing
end

function simd_mult!(w, a, X, y)
    @assert length(w) == size(X, 1)
    @assert length(y) == size(X, 2)
    @simd for i=1:length(w)
        w[i] = 0
        @simd for j=1:length(y)
            @inbounds w[i] += a * X[i,j] * y[j]
        end
    end
    nothing
end
