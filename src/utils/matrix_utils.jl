
function istriustrict(A::Matrix)
    m, n = size(A)
    @inbounds for j = 1:min(n,m-1), i = j:m
        if A[i,j] != 0.
            return false
        end
    end
    return true
end

function istrilstrict(A::Matrix)
    m, n = size(A)
    @inbounds for j = 2:n, i = 1:min(j,m)
        if A[i,j] != 0.
            return false
        end
    end
    return true
end

function simd_scale!(x, a)
    @inbounds for i=1:length(x)
        x[i] *= a
    end
    nothing
end

function simd_copy!(x, y)
    @assert length(x) == length(y)
    @inbounds for i=1:length(y)
        y[i] = x[i]
    end
    nothing
end

function simd_copy_scale!(a, x, y)
    @assert length(x) == length(y)
    @inbounds for i=1:length(y)
        y[i] = a*x[i]
    end
    nothing
end

"Copy the first dimension of a 2D array y into a 1D array x."
function simd_copy_xy_first!(x, y, j)
    @assert length(x) == size(y, 1)
    @inbounds for i=1:length(x)
        x[i] = y[i,j]
    end
    nothing
end

"Copy the first dimension of a 3D array y into a 1D array x."
function simd_copy_xy_first!(x, y, j, k)
    @assert length(x) == size(y, 1)
    @inbounds for i=1:length(x)
        x[i] = y[i,j,k]
    end
    nothing
end

"Copy a 1D array x into the first dimension of a 2D array y."
function simd_copy_yx_first!(x, y, j)
    @assert length(x) == size(y, 1)
    @inbounds for i=1:size(y, 1)
        y[i,j] = x[i]
    end
    nothing
end

"Copy a 1D array x into the first dimension of a 3D array y."
function simd_copy_yx_first!(x, y, j, k)
    @assert length(x) == size(y, 1)
    @inbounds for i=1:size(y, 1)
        y[i,j,k] = x[i]
    end
    nothing
end

function simd_copy_yx_first_last!(x, y, j)
    @assert size(x, 1) == size(y, 1)
    @assert size(x, 2) == size(y, 3)
    @inbounds for k=1:size(y, 3)
        for i=1:size(y, 1)
            y[i,j,k] = x[i,k]
        end
    end
    nothing
end

function simd_xpy!(x, y)
    @assert length(x) == length(y)
    @inbounds for i=1:length(y)
        y[i] += x[i]
    end
    nothing
end

function simd_axpy!(a, x, y)
    @assert length(x) == length(y)
    @inbounds for i=1:length(y)
        y[i] += a*x[i]
    end
    nothing
end

function simd_wxpy!(w, x, y)
    @assert length(x) == length(y) == length(w)
    @inbounds for i=1:length(w)
        w[i] = x[i] + y[i]
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
    @inbounds for i=1:length(w)
        w[i] = 0.
        for j=1:length(y)
            w[i] += X[i,j] * y[j]
        end
    end
    nothing
end

function simd_mult!(w, a, X, y)
    @assert length(w) == size(X, 1)
    @assert length(y) == size(X, 2)
    @inbounds for i=1:length(w)
        w[i] = 0
        for j=1:length(y)
            w[i] += a * X[i,j] * y[j]
        end
    end
    nothing
end
