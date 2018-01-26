
abstract type DataSeries{T,N} <: AbstractArray{T,N} end


struct SDataSeries{T,N} <: DataSeries{T,N}
    nd::Int
    nt::Int
    ni::Int
    d::Array{T,N}

    function SDataSeries{T,N}(nd, nt, ni) where {T,N}
        @assert T <: Number
        @assert nd > 0
        @assert nt ≥ 0
        @assert ni > 0

        @assert N ∈ (1,2,3)

        if N == 1
            d = zeros(T, nt+1)
        elseif N == 2
            d = zeros(T, nd, nt+1)
        elseif N == 3
            d = zeros(T, nd, nt+1, ni)
        end

        new(nd, nt, ni, d)
    end

    function SDataSeries{T,N}(nd, nt, ni, d) where {T,N}
        new(nd, nt, ni, d)
    end
end

function SDataSeries(T, nd::Int, nt::Int, ni::Int)
    # ni == 1 ? (nd == 1 ? N = 1 : N = 2) : N = 3
    ni == 1 ? N = 2 : N = 3
    return SDataSeries{T,N}(nd, nt, ni)
end

function SDataSeries(d::Array{T,1}) where {T}
    nd = 1
    nt = size(d,1)
    ni = 1
    return SDataSeries{T,1}(nd, nt, ni, d)
end

function SDataSeries(d::Array{T,2}) where {T}
    nd = size(d,1)
    nt = size(d,2)
    ni = 1
    return SDataSeries{T,2}(nd, nt, ni, d)
end

function SDataSeries(d::Array{T,3}) where {T}
    nd = size(d,1)
    nt = size(d,2)
    ni = size(d,3)
    return SDataSeries{T,3}(nd, nt, ni, d)
end


struct PDataSeries{T,N} <: DataSeries{T,N}
    nd::Int
    nt::Int
    ni::Int
    d::SharedArray{T,N}

    function PDataSeries{T,N}(nd, nt, ni) where {T,N}
        @assert T <: Number
        @assert nd > 0
        @assert nt ≥ 0
        @assert ni > 0

        @assert N ∈ (1,2,3)

        if N == 1
            d = SharedArray(T, (nt+1))
        elseif N == 2
            d = SharedArray(T, (nd, nt+1))
        elseif N == 3
            d = SharedArray(T, (nd, nt+1, ni))
        end

        new(nd, nt, ni, d)
    end

    function PDataSeries{T,N}(nd, nt, ni, d::SharedArray{T,N}) where {T,N}
        new(nd, nt, ni, d)
    end

    function PDataSeries{T,N}(nd, nt, ni, d::Array{T,N}) where {T,N}
        sd = SharedArray(eltype(d),size(d));
        copy!(sd,d)
        new(nd, nt, ni, sd)
    end
end

function PDataSeries(T, nd::Int, nt::Int, ni::Int)
    # ni == 1 ? (nd == 1 ? N = 1 : N = 2) : N = 3
    ni == 1 ? N = 2 : N = 3
    return PDataSeries{T,N}(nd, nt, ni)
end

function PDataSeries(d::Array{T,1}) where {T}
    nd = 1
    nt = size(d,1)
    ni = 1
    return PDataSeries{T,1}(nd, nt, ni, d)
end

function PDataSeries(d::Array{T,2}) where {T}
    nd = size(d,1)
    nt = size(d,2)
    ni = 1
    return PDataSeries{T,2}(nd, nt, ni, d)
end

function PDataSeries(d::Array{T,3}) where {T}
    nd = size(d,1)
    nt = size(d,2)
    ni = size(d,3)
    return PDataSeries{T,3}(nd, nt, ni, d)
end


function similar(ds::DataSeries{T,N}) where {T,N}
    typeof(ds){T,N}(ds.nd, ds.nt, ds.ni)
end

Base.eltype(ds::DataSeries{T,N}) where {T,N} = T
Base.ndims(ds::DataSeries{T,N}) where {T,N} = N
Base.size(ds::DataSeries) = size(ds.d)
Base.length(ds::DataSeries) = length(ds.d)
Base.endof(ds::DataSeries{T,1}) where {T} = (ds.nt)
Base.endof(ds::DataSeries{T,2}) where {T} = (ds.nd, ds.nt)
Base.endof(ds::DataSeries{T,3}) where {T} = (ds.nd, ds.nt, ds.ni)
Base.indices(ds::DataSeries{T,1}) where {T} = (0:ds.nt)
Base.indices(ds::DataSeries{T,2}) where {T} = (1:ds.nd, 0:ds.nt)
Base.indices(ds::DataSeries{T,3}) where {T} = (1:ds.nd, 0:ds.nt, 1:ds.ni)
Base.strides(ds::DataSeries) = (strides(ds.d))

function get_data!(ds::DataSeries{T,1}, n, k=1) where {T}
    j = n+1
    @assert length(x) == size(ds.d, 1)
    @assert j ≤ size(ds.d, 2)
    @assert k == 1
    @inbounds return ds.d[j]
end

function get_data!(ds::DataSeries{T,2}, x::Union{Array{T,1}, Array{Double{T},1}}, n, k=1) where {T}
    j = n+1
    @assert j ≤ size(ds.d, 2)
    @assert k == 1
    @inbounds for i in eachindex(x)
        x[i] = ds.d[i,j]
    end
end

function get_data!(ds::DataSeries{T,3}, x::Union{Array{T,2}, Array{Double{T},2}}, n) where {T}
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

function get_data!(ds::DataSeries{T,3}, x::Union{Array{T,1}, Array{Double{T},1}}, n, k) where {T}
    j = n+1
    @assert size(x,1) == size(ds.d, 1)
    @assert j ≤ size(ds.d, 2)
    @assert k ≤ size(ds.d, 3)
    @inbounds for i in eachindex(x)
        x[i] = ds.d[i,j,k]
    end
end

function set_data!(ds::DataSeries{T,1}, x::Union{T, Double{T}}, n, k=1) where {T}
    j = n+1
    @assert j ≤ size(ds.d, 2)
    @assert k == 1
    @inbounds ds.d[j] = x[i]
end

function set_data!(ds::DataSeries{T,2}, x::Union{Array{T,1}, Array{Double{T},1}}, n, k=1) where {T}
    j = n+1
    @assert length(x) == size(ds.d, 1)
    @assert j ≤ size(ds.d, 2)
    @assert k == 1
    @inbounds for i in 1:size(ds.d, 1)
        ds.d[i,j] = x[i]
    end
end

function set_data!(ds::DataSeries{T,3}, x::Union{Array{T,2}, Array{Double{T},2}}, n) where {T}
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

function set_data!(ds::DataSeries{T,3}, x::Union{Array{T,1}, Array{Double{T},1}}, n, k) where {T}
    j = n+1
    @assert size(x,1) == size(ds.d, 1)
    @assert j ≤ size(ds.d, 2)
    @assert k ≤ size(ds.d, 3)
    @inbounds for i in 1:size(ds.d, 1)
        ds.d[i,j,k] = x[i]
    end
end

function reset!(ds::DataSeries{T,1}) where {T}
    @inbounds ds[0] = ds[end]
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

@inline function Base.getindex(ds::DataSeries{T,1}, j::Int) where {T}
    @boundscheck checkbounds(ds.d, j+1)
    @inbounds r = getindex(ds.d, j+1)
    return r
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

@inline function Base.setindex!(ds::DataSeries{T,1}, x, j::Int) where {T}
    @assert length(x) == 1
    @boundscheck checkbounds(ds.d, j+1)
    @inbounds setindex!(ds.d, x, j+1)
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
