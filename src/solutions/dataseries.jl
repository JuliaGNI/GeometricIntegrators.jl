
using Base: Indices, tail

abstract type DataSeries{T,N} <: AbstractArray{T,N} end

for (TDataSeries, TArray) in
    ((:SDataSeries, :Array),
     (:PDataSeries, :SharedArray))
    @eval begin
        struct $TDataSeries{T,N} <: DataSeries{T,N}
            nd::Int
            nt::Int
            ni::Int
            d::$TArray{T,N}

            function $TDataSeries{T,N}(nd, nt, ni) where {T,N}
                @assert T <: Number
                @assert nd > 0
                @assert nt ≥ 0
                @assert ni > 0

                @assert N ∈ (1,2,3)

                if N == 1
                    d = $TArray{T}(undef, nt+1)
                elseif N == 2
                    d = $TArray{T}(undef, nd, nt+1)
                elseif N == 3
                    d = $TArray{T}(undef, nd, nt+1, ni)
                end

                fill!(d,zero(T))

                new(nd, nt, ni, d)
            end
        end


        function $TDataSeries(T, nt::Int)
            return $TDataSeries{T,1}(1, nt, 1)
        end

        function $TDataSeries(T, nd::Int, nt::Int)
            return $TDataSeries{T,2}(nd, nt, 1)
        end

        function $TDataSeries(T, nd::Int, nt::Int, ni::Int)
            ni == 1 ? N = 2 : N = 3
            return $TDataSeries{T,N}(nd, nt, ni)
        end

        function $TDataSeries(d::AbstractArray{T,1}) where {T}
            nd = 1
            nt = size(d,1)-1
            ni = 1
            ds = $TDataSeries{T,1}(nd, nt, ni)
            copy!(ds.d, d[1,:,1])
            return ds
        end

        function $TDataSeries(d::AbstractArray{T,2}) where {T}
            nd = size(d,1)
            nt = size(d,2)-1
            ni = 1
            ds = $TDataSeries{T,2}(nd, nt, ni)
            copy!(ds.d, d[:,:,1])
            return ds
        end

        function $TDataSeries(d::AbstractArray{T,3}) where {T}
            nd = size(d,1)
            nt = size(d,2)-1
            ni = size(d,3)
            ds = $TDataSeries{T,3}(nd, nt, ni)
            copy!(ds.d, d)
            return ds
        end
    end
end

Base.:(==)(ds1::DataSeries, ds2::DataSeries) = (
                                ds1.nd == ds2.nd
                             && ds1.nt == ds2.nt
                             && ds1.ni == ds2.ni
                             && ds1.d  == ds2.d)

function Base.show(io::IO, ds::DataSeries{T,N}) where {T,N}
    print(io, "DataSeries with data type ", T, " and ", N, " dimensions:\n")
    print(io, ds.d)
end

function similar(ds::DataSeries{T,N}) where {T,N}
    typeof(ds){T,N}(ds.nd, ds.nt, ds.ni)
end

Base.parent(ds::DataSeries) = ds.d
Base.eltype(ds::DataSeries{T,N}) where {T,N} = T
Base.ndims(ds::DataSeries{T,N}) where {T,N} = N

# errmsg_size(ds::DataSeries) = error("size not supported for data series with axes $(axes(ds))")
# Base.size(ds::DataSeries) = errmsg_size(ds)
# Base.size(ds::DataSeries, d) = errmsg_size(ds)
Base.size(ds::DataSeries) = size(ds.d)
Base.size(ds::DataSeries, d) = size(ds.d, d)

Base.eachindex(::IndexCartesian, ds::DataSeries) = CartesianIndices(axes(ds))
Base.eachindex(::IndexLinear, ds::DataSeries) = axes(ds, 1)

Base.firstindex(ds::DataSeries{T,1}) where {T}   = 0
Base.firstindex(ds::DataSeries{T,N}) where {T,N} = 1

Base.firstindex(ds::DataSeries{T,1}, d) where {T} = d ≥ 1 && d ≤ 1 ? 0 : 1
Base.firstindex(ds::DataSeries{T,2}, d) where {T} = d ≥ 1 && d ≤ 2 ? (1,0)[d] : 1
Base.firstindex(ds::DataSeries{T,3}, d) where {T} = d ≥ 1 && d ≤ 3 ? (1,0,1)[d] : 1

Base.lastindex(ds::DataSeries{T,1}) where {T} = ds.nt
Base.lastindex(ds::DataSeries{T,2}) where {T} = ds.nd*(ds.nt+1)
Base.lastindex(ds::DataSeries{T,3}) where {T} = ds.nd*(ds.nt+1)*ds.ni

Base.lastindex(ds::DataSeries{T,1}, d) where {T} = d ≥ 1 && d ≤ 1 ? ds.nt : 1
Base.lastindex(ds::DataSeries{T,2}, d) where {T} = d ≥ 1 && d ≤ 2 ? (ds.nd, ds.nt)[d] : 1
Base.lastindex(ds::DataSeries{T,3}, d) where {T} = d ≥ 1 && d ≤ 3 ? (ds.nd, ds.nt, ds.ni)[d] : 1

@inline Base.axes(ds::DataSeries{T,1}) where {T} = (0:ds.nt,)
@inline Base.axes(ds::DataSeries{T,2}) where {T} = (1:ds.nd, 0:ds.nt)
@inline Base.axes(ds::DataSeries{T,3}) where {T} = (1:ds.nd, 0:ds.nt, 1:ds.ni)
@inline Base.axes(ds::DataSeries{T,N}, d) where {T,N} = d ≥ 1 && d ≤ N ? axes(ds)[d] : (1:1)

Base.strides(ds::DataSeries) = strides(ds.d)


function get_data!(ds::DataSeries{T,1}, n, k=1) where {T}
    @assert n ≥ 0 && n ≤ ds.nt
    @assert k == 1
    @inbounds return ds.d[n+1]
end

function get_data!(ds::DataSeries{T,2}, x::Union{Array{T,1}, Array{TwicePrecision{T},1}}, n, k=1) where {T}
    @assert size(x,1) == ds.nd
    @assert n ≥ 0 && n ≤ ds.nt
    @assert k == 1
    j = n+1
    @inbounds for i in eachindex(x)
        x[i] = ds.d[i,j]
    end
end

function get_data!(ds::DataSeries{T,3}, x::Union{Array{T,2}, Array{TwicePrecision{T},2}}, n) where {T}
    @assert size(x,1) == ds.nd
    @assert size(x,2) == ds.ni
    @assert n ≥ 0 && n ≤ ds.nt
    j = n+1
    @inbounds for k in axes(x, 2)
        for i in axes(x, 1)
            x[i,k] = ds.d[i,j,k]
        end
    end
end

function get_data!(ds::DataSeries{T,3}, x::Union{Array{T,1}, Array{TwicePrecision{T},1}}, n, k) where {T}
    @assert size(x,1) == ds.nd
    @assert n ≥ 0 && n ≤ ds.nt
    @assert k ≥ 1 && k ≤ ds.ni
    j = n+1
    @inbounds for i in eachindex(x)
        x[i] = ds.d[i,j,k]
    end
end

function set_data!(ds::DataSeries{T,1}, x::Union{T, TwicePrecision{T}}, n, k=1) where {T}
    @assert n ≥ 0 && n ≤ ds.nt
    @assert k == 1
    j = n+1
    @inbounds ds.d[j] = x[i]
end

function set_data!(ds::DataSeries{T,2}, x::Union{Array{T,1}, Array{TwicePrecision{T},1}}, n, k=1) where {T}
    @assert size(x,1) == ds.nd
    @assert n ≥ 0 && n ≤ ds.nt
    @assert k == 1
    j = n+1
    @inbounds for i in axes(ds.d, 1)
        ds.d[i,j] = x[i]
    end
end

function set_data!(ds::DataSeries{T,3}, x::Union{Array{T,2}, Array{TwicePrecision{T},2}}, n) where {T}
    @assert size(x,1) == ds.nd
    @assert size(x,2) == ds.ni
    @assert n ≥ 0 && n ≤ ds.nt
    j = n+1
    @inbounds for k in axes(ds.d, 3)
        for i in axes(ds.d, 1)
            ds.d[i,j,k] = x[i,k]
        end
    end
end

function set_data!(ds::DataSeries{T,3}, x::Union{Array{T,1}, Array{TwicePrecision{T},1}}, n, k) where {T}
    @assert size(x,1) == ds.nd
    @assert n ≥ 0 && n ≤ ds.nt
    @assert k ≥ 1 && k ≤ ds.ni
    j = n+1
    @inbounds for i in axes(ds.d, 1)
        ds.d[i,j,k] = x[i]
    end
end

function CommonFunctions.reset!(ds::DataSeries{T,1}) where {T}
    @inbounds ds[0] = ds[end]
end

function CommonFunctions.reset!(ds::DataSeries{T,2}) where {T}
    @inbounds for i in axes(ds, 1)
        ds[i,0] = ds[i,end]
    end
end

function CommonFunctions.reset!(ds::DataSeries{T,3}) where {T}
    @inbounds for k in axes(ds.d, 3)
        for i in axes(ds.d, 1)
            ds[i,0,k] = ds[i,end,k]
        end
    end
end

@inline function Base.getindex(ds::DataSeries{T,1}, j::Union{UnitRange,Int}) where {T}
    @boundscheck checkbounds(ds.d, j.+1)
    @inbounds r = getindex(ds.d, j.+1)
    return r
end

@inline function Base.getindex(ds::DataSeries{T,2}, i::Union{UnitRange,Int}, j::Union{UnitRange,Int}) where {T}
    @boundscheck checkbounds(ds.d, i, j.+1)
    @inbounds r = getindex(ds.d, i, j.+1)
    return r
end

@inline function Base.getindex(ds::DataSeries{T,2}, j::Union{UnitRange,Int}) where {T}
    @boundscheck checkbounds(ds.d, :, j.+1)
    @inbounds r = getindex(ds.d, :, j+1)
    return r
end

@inline function Base.getindex(ds::DataSeries{T,3}, i::Union{UnitRange,Int}, j::Union{UnitRange,Int}, k::Union{UnitRange,Int}) where {T}
    @boundscheck checkbounds(ds.d, i, j.+1, k)
    @inbounds r = getindex(ds.d, i, j.+1, k)
    return r
end

@inline function Base.getindex(ds::DataSeries{T,3}, j::Union{UnitRange,Int}, k::Union{UnitRange,Int}) where {T}
    @boundscheck checkbounds(ds.d, :, j.+1, k)
    @inbounds r = getindex(ds.d, :, j.+1, k)
    return r
end

@inline function Base.getindex(ds::DataSeries{T,3}, k::Union{UnitRange,Int}) where {T}
    @boundscheck checkbounds(ds.d, :, :, k)
    @inbounds r = getindex(ds.d, :, :, k)
    return r
end

@inline function Base.setindex!(ds::DataSeries{T,1}, x, j::Union{UnitRange,Int}) where {T}
    @assert length(x) == 1
    @boundscheck checkbounds(ds.d, j.+1)
    @inbounds setindex!(ds.d, x, j.+1)
end

@inline function Base.setindex!(ds::DataSeries{T,2}, x, i::Union{UnitRange,Int}, j::Union{UnitRange,Int}) where {T}
    @assert length(x) == 1
    @boundscheck checkbounds(ds.d, i, j.+1)
    @inbounds setindex!(ds.d, x, i, j.+1)
end

@inline function Base.setindex!(ds::DataSeries{T,2}, x, j::Union{UnitRange,Int}) where {T}
    @assert ndims(x) == 1
    @assert length(x) == size(ds.d, 1)
    @boundscheck checkbounds(ds.d, :, j.+1)
    @inbounds setindex!(ds.d, x, :, j.+1)
end

@inline function Base.setindex!(ds::DataSeries{T,3}, x, i::Union{UnitRange,Int}, j::Union{UnitRange,Int}, k::Union{UnitRange,Int}) where {T}
    @assert length(x) == 1
    @boundscheck checkbounds(ds.d, i, j.+1, k)
    @inbounds setindex!(ds.d, x, i, j.+1, k)
end

@inline function Base.setindex!(ds::DataSeries{T,3}, x, j::Union{UnitRange,Int}, k::Union{UnitRange,Int}) where {T}
    @assert ndims(x) == 1
    @assert length(x) == size(ds.d, 1)
    @boundscheck checkbounds(ds.d, :, j.+1, k)
    @inbounds setindex!(ds.d, x, :, j.+1, k)
end

@inline function Base.setindex!(ds::DataSeries{T,3}, x, k::Union{UnitRange,Int}) where {T}
    @assert ndims(x) == 2
    @assert size(x,1) == size(ds.d, 1)
    @assert size(x,2) == size(ds.d, 2)
    @boundscheck checkbounds(ds.d, :, :, k)
    @inbounds setindex!(ds.d, x, :, :, k)
end

# TODO Implement convert() to/from array.
