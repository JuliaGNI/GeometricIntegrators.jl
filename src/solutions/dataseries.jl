
abstract type DataSeries{T,N} <: AbstractArray{T,N} end


struct SDataSeries{T,N} <: DataSeries{T,N}
    nd::Int
    nt::Int
    ni::Int
    d::Array{T,N}

    function SDataSeries{T,N}(nd, nt, ni) where {T,N}
        @assert T <: Number
        @assert nd > 0
        @assert nt > 0
        @assert ni > 0

        @assert N ∈ (2,3)

        if N == 2
            d = zeros(T, nd, nt+1)
        elseif N == 3
            d = zeros(T, nd, nt+1, ni)
        end

        new(nd, nt, ni, d)
    end
end

function SDataSeries(T, nd::Int, nt::Int, ni::Int)
    ni == 1 ? N = 2 : N = 3
    return SDataSeries{T,N}(nd, nt, ni)
end



struct PDataSeries{T,N} <: DataSeries{T,N}
    nd::Int
    nt::Int
    ni::Int
    d::SharedArray{T,N}

    function PDataSeries{T,N}(nd, nt, ni) where {T,N}
        @assert T <: Number
        @assert nd > 0
        @assert nt > 0
        @assert ni > 0

        @assert N ∈ (2,3)

        if N == 2
            d = SharedArray(T, (nd, nt+1))
        elseif N == 3
            d = SharedArray(T, (nd, nt+1, ni))
        end

        new(nd, nt, ni, d)
    end
end

function PDataSeries(T, nd::Int, nt::Int, ni::Int)
    ni == 1 ? N = 2 : N = 3
    return PDataSeries{T,N}(nd, nt, ni)
end


function similar(ds::DataSeries{T,N}) where {T,N}
    typeof(ds){T,N}(ds.nd, ds.nt, ds.ni)
end

Base.eltype(ds::DataSeries{T,N}) where {T,N} = T
Base.ndims(ds::DataSeries{T,N}) where {T,N} = N
Base.size(ds::DataSeries) = size(ds.d)
Base.length(ds::DataSeries) = length(ds.d)
Base.endof(ds::DataSeries{T,2}) where {T} = (ds.nd, ds.nt)
Base.endof(ds::DataSeries{T,3}) where {T} = (ds.nd, ds.nt, ds.ni)
Base.indices(ds::DataSeries{T,2}) where {T} = (1:ds.nd, 0:ds.nt)
Base.indices(ds::DataSeries{T,3}) where {T} = (1:ds.nd, 0:ds.nt, 1:ds.ni)
Base.strides(ds::DataSeries) = (strides(ds.d))

function get_data!(ds::DataSeries{T,2}, x::Array{T,1}, n, k=1) where {T}
    j = n+1
    @assert length(x) == size(ds.d, 1)
    @assert j ≤ size(ds.d, 2)
    @assert k == 1
    @inbounds for i in eachindex(x)
        x[i] = ds.d[i,j]
    end
end

function get_data!(ds::DataSeries{T,3}, x::Array{T,2}, n) where {T}
    j = n+1
    @assert size(x,1) == size(ds.d, 1)
    @assert size(x,2) == size(ds.d, 3)
    @assert j ≤ size(ds.d, 2)
    @inbounds for k in 1:size(x, 2)
        for i in 1:size(x, 1)
            x[i,k] = ds.d[i,j,k]
        end
    end
end

function get_data!(ds::DataSeries{T,3}, x::Array{T,1}, n, k) where {T}
    j = n+1
    @assert size(x,1) == size(ds.d, 1)
    @assert j ≤ size(ds.d, 2)
    @assert k ≤ size(ds.d, 3)
    @inbounds for i in eachindex(x)
        x[i] = ds.d[i,j,k]
    end
end

function set_data!(ds::DataSeries{T,2}, x::Array{T,1}, n, k=1) where {T}
    j = n+1
    @assert length(x) == size(ds.d, 1)
    @assert j ≤ size(ds.d, 2)
    @assert k == 1
    @inbounds for i in 1:size(ds.d, 1)
        ds.d[i,j] = x[i]
    end
end

function set_data!(ds::DataSeries{T,3}, x::Array{T,2}, n) where {T}
    j = n+1
    @assert size(x,1) == size(ds.d, 1)
    @assert size(x,2) == size(ds.d, 3)
    @assert j ≤ size(ds.d, 2)
    @inbounds for k in 1:size(ds.d, 3)
        for i in 1:size(ds.d, 1)
            ds.d[i,j,k] = x[i,k]
        end
    end
end

function set_data!(ds::DataSeries{T,3}, x::Array{T,1}, n, k) where {T}
    j = n+1
    @assert size(x,1) == size(ds.d, 1)
    @assert j ≤ size(ds.d, 2)
    @assert k ≤ size(ds.d, 3)
    @inbounds for i in 1:size(ds.d, 1)
        ds.d[i,j,k] = x[i]
    end
end

function reset!(ds::DataSeries{T,2}) where {T}
    @inbounds for i in 1:size(ds,1)
        ds[i,0] = ds[i,end]
    end
end

function reset!(ds::DataSeries{T,3}) where {T}
    @inbounds for k in 1:size(ds,3)
        for i in 1:size(ds,1)
            ds[i,0,k] = ds[i,end,k]
        end
    end
end

@inline function Base.getindex(ds::DataSeries{T,2}, i::Int, j::Int) where {T}
    @boundscheck checkbounds(ds.d, i, j+1)
    @inbounds r = getindex(ds.d, i, j+1)
    return r
end

@inline function Base.getindex(ds::DataSeries{T,2}, j::Int) where {T}
    @boundscheck checkbounds(ds.d, :, j+1)
    @inbounds r = getindex(ds.d, :, j+1)
    return r
end

@inline function Base.getindex(ds::DataSeries{T,3}, i::Int, j::Int, k::Int) where {T}
    @boundscheck checkbounds(ds.d, i, j+1, k)
    @inbounds r = getindex(ds.d, i, j+1, k)
    return r
end

@inline function Base.getindex(ds::DataSeries{T,3}, j::Int, k::Int) where {T}
    @boundscheck checkbounds(ds.d, :, j+1, k)
    @inbounds r = getindex(ds.d, :, j+1, k)
    return r
end

@inline function Base.getindex(ds::DataSeries{T,3}, k::Int) where {T}
    @boundscheck checkbounds(ds.d, :, :, k)
    @inbounds r = getindex(ds.d, :, :, k)
    return r
end

@inline function Base.setindex!(ds::DataSeries{T,2}, x, i::Int, j::Int) where {T}
    @assert length(x) == 1
    @boundscheck checkbounds(ds.d, i, j+1)
    @inbounds setindex!(ds.d, x, i, j+1)
end

@inline function Base.setindex!(ds::DataSeries{T,2}, x, j::Int) where {T}
    @assert ndims(x) == 1
    @assert length(x) == size(ds.d, 1)
    @boundscheck checkbounds(ds.d, :, j+1)
    @inbounds setindex!(ds.d, x, :, j+1)
end

@inline function Base.setindex!(ds::DataSeries{T,3}, x, i::Int, j::Int, k::Int) where {T}
    @assert length(x) == 1
    @boundscheck checkbounds(ds.d, i, j+1, k)
    @inbounds setindex!(ds.d, x, i, j+1, k)
end

@inline function Base.setindex!(ds::DataSeries{T,3}, x, j::Int, k::Int) where {T}
    @assert ndims(x) == 1
    @assert length(x) == size(ds.d, 1)
    @boundscheck checkbounds(ds.d, :, j+1, k)
    @inbounds setindex!(ds.d, x, :, j+1, k)
end

@inline function Base.setindex!(ds::DataSeries{T,3}, x, k::Int) where {T}
    @assert ndims(x) == 2
    @assert size(x,1) == size(ds.d, 1)
    @assert size(x,2) == size(ds.d, 2)
    @boundscheck checkbounds(ds.d, :, :, k)
    @inbounds setindex!(ds.d, x, :, :, k)
end

# TODO Implement convert() to/from array.
