
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


"Copy the first dimension of a 2D array y into a 1D array x."
function simd_copy_xy_first!(x::Array{T,1}, y::Array{T,2}, j) where {T}
    @assert length(x) == size(y, 1)
    @inbounds for i in eachindex(x)
        x[i] = y[i,j]
    end
    nothing
end

"Copy the first dimension of a 3D array y into a 1D array x."
function simd_copy_xy_first!(x::Array{T,1}, y::Array{T,3}, j, k) where {T}
    @assert length(x) == size(y, 1)
    @inbounds for i in eachindex(x)
        x[i] = y[i,j,k]
    end
    nothing
end

"Copy the first dimension of a 3D array y into a 1D array x."
function simd_copy_xy_first!(x::Array{T,2}, y::Array{T,3}, k) where {T}
    @assert size(x,1) == size(y, 1)
    @assert size(x,2) == size(y, 2)
    @inbounds for i in 1:size(x,1)
        for j=1:size(x,2)
            x[i,j] = y[i,j,k]
        end
    end
    nothing
end

"Copy a 1D array x into the first dimension of a 2D array y."
function simd_copy_yx_first!(x::Array{T,1}, y::Array{T,2}, j) where {T}
    @assert length(x) == size(y, 1)
    @inbounds for i in 1:size(y, 1)
        y[i,j] = x[i]
    end
    nothing
end

"Copy a 1D array x into the second dimension of a 2D array y."
function simd_copy_yx_second!(x::Array{T,1}, y::Array{T,2}, j) where {T}
    @assert length(x) == size(y, 2)
    @inbounds for i in 1:size(y, 2)
        y[j,i] = x[i]
    end
    nothing
end

"Copy a 1D array x into the first dimension of a 3D array y."
function simd_copy_yx_first!(x::Array{T,1}, y::Array{T,3}, j, k) where {T}
    @assert length(x) == size(y, 1)
    @inbounds for i in 1:size(y, 1)
        y[i,j,k] = x[i]
    end
    nothing
end

"Copy a 1D array x into the first dimension of a 3D array y."
function simd_copy_yx_first!(x::Array{T,2}, y::Array{T,3}, k) where {T}
    @assert size(x, 1) == size(y, 1)
    @assert size(x, 2) == size(y, 2)
    @inbounds for j in 1:size(y, 2)
        for i in 1:size(y, 1)
            y[i,j,k] = x[i,j]
        end
    end
    nothing
end

function simd_copy_yx_first_last!(x::Array{T,2}, y::Array{T,3}, j) where {T}
    @assert size(x, 1) == size(y, 1)
    @assert size(x, 2) == size(y, 3)
    @inbounds for k=1:size(y, 3)
        for i in 1:size(y, 1)
            y[i,j,k] = x[i,k]
        end
    end
    nothing
end

function simd_axpy!(a::T, x::DenseMatrix{T}, y::DenseMatrix{T}, e::DenseMatrix{T}) where {T}
    @assert size(x) == size(y) == size(e)

    local err::T
    local ty::T

    @inbounds for i in eachindex(x,y,e)
        err = e[i] + a*x[i]
        ty  = y[i]
        y[i] = err + ty
        e[i] = err + (ty - y[i])
    end
    nothing
end


function simd_aXbpy!(a, b::Vector{T}, X::Matrix{T}, y::Vector{T}) where {T}
    @assert length(y) == size(X, 1)
    @assert length(b) == size(X, 2)
    local ty::T
    @inbounds for i in eachindex(y)
        ty = zero(T)
        for j=eachindex(b)
            ty += X[i,j] * b[j]
        end
        y[i] = a*ty
    end
    nothing
end

function simd_abXpy!(a::T, b::Vector{T}, X::Matrix{T}, y::Vector{T}) where {T}
    @assert length(y) == size(X, 2)
    @assert length(b) == size(X, 1)
    local ty::T
    @inbounds for i in eachindex(y)
        ty = zero(T)
        for j=eachindex(b)
            ty += b[j] * X[j,i]
        end
        y[i] += a*ty
    end
    nothing
end

function simd_mult!(w::Vector{T}, X::Matrix{T}, y::Vector{T}) where {T}
    @assert length(w) == size(X, 1)
    @assert length(y) == size(X, 2)
    local tw::T
    @inbounds for i in eachindex(w)
        tw = zero(T)
        for j=eachindex(y)
            tw += X[i,j] * y[j]
        end
        w[i] = tw
    end
    nothing
end

function simd_mult!(w::Vector{T}, y::Vector{T}, X::Matrix{T}) where {T}
    @assert length(w) == size(X, 2)
    @assert length(y) == size(X, 1)
    local tw::T
    @inbounds for i in eachindex(w)
        tw = 0
        for j=eachindex(y)
            tw += y[j] * X[j,i]
        end
        w[i] = tw
    end
    nothing
end

function simd_mult!(W::Matrix{T}, X::Matrix{T}, Y::Matrix{T}) where {T}
    @assert size(W,1) == size(X, 1)
    @assert size(W,2) == size(Y, 2)
    @assert size(X,2) == size(Y, 1)
    local tw::T
    @inbounds for j in 1:size(W,2)
        for i in 1:size(W,1)
            tw = 0
            for k in 1:size(X,2)
                tw += X[i,k] * Y[k,j]
            end
            W[i,j] = tw
        end
    end
    nothing
end
