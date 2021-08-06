
using Base: Indices

abstract type DataSeries{T <: Union{Number,AbstractArray}, N} <: AbstractArray{T,N} end

function initialize!(ds::DataSeries{T,1}, q₀::T) where {T}
    for j in eachindex(ds)
        ds[j] = zero(q₀)
    end
    return ds
end

function initialize!(ds::DataSeries{T,1}, q₀::Vector{T}) where {T}
    @assert nsamples(ds) == length(q₀) == 1
    initialize!(ds, q₀[begin])
end

function initialize!(ds::DataSeries{T,2}, q₀::Vector{T}) where {T}
    @assert nsamples(ds) == length(q₀)
    for k in axes(ds,2)
        for j in axes(ds,1)
            ds[j,k] = zero(q₀[k])
        end
    end
    return ds
end


for (TDataSeries, TArray) in
    ((:SDataSeries, :Array),
     (:PDataSeries, :SharedArray))
    @eval begin
        struct $TDataSeries{T <: Union{Number,AbstractArray},N} <: DataSeries{T,N}
            d::$TArray{T,N}

            function $TDataSeries{T,N}(nt, ni) where {T <: Union{Number,AbstractArray}, N}
                @assert nt ≥ 0
                @assert ni > 0

                @assert N ∈ (1,2)

                if $TArray == Array
                    if N == 1
                        d = $TArray{T,N}(undef, nt+1)
                    elseif N == 2
                        d = $TArray{T,N}(undef, nt+1, ni)
                    end
                elseif $TArray == SharedArray
                    if N == 1
                        d = $TArray{T,N}(nt+1)
                    elseif N == 2
                        d = $TArray{T,N}(nt+1, ni)
                    end
                end

                new(d)
            end
        end

        function $TDataSeries(::Type{T}, nt::Int, ni::Int=1) where {T <: Number}
            ds = $TDataSeries{T, ni == 1 ? 1 : 2}(nt, ni)
            initialize!(ds, zero(T))
        end

        function $TDataSeries(::Type{T}, nd::Int, nt::Int, ni::Int) where {T <: Number}
            ds = $TDataSeries{Vector{T}, ni == 1 ? 1 : 2}(nt, ni)
            initialize!(ds, [zeros(T, nd) for _ in 1:ni])
        end

        function $TDataSeries(q₀::T, nt::Int, ni::Int=1) where {T <: Union{Number,AbstractArray{<:Number}}}
            @assert ni == 1
            ds = $TDataSeries{T,1}(nt, ni)
            initialize!(ds, q₀)
        end

        function $TDataSeries(q₀::Vector{T}, nt::Int, ni::Int=1) where {T <: Union{Number}}
            @assert ni == 1 || ni == length(q₀)
            ds = ( ni == 1 ? $TDataSeries{Vector{T},1}(nt, ni) : $TDataSeries{T,2}(nt, ni) )
            initialize!(ds, q₀)
        end

        function $TDataSeries(q₀::Vector{T}, nt::Int, ni::Int=length(q₀)) where {T <: AbstractArray{<:Number}}
            @assert ni == length(q₀)
            ds = $TDataSeries{T, ni == 1 ? 1 : 2}(nt, ni)
            initialize!(ds, q₀)
        end

        function $TDataSeries(d::AbstractArray{T,1}) where {T}
            nt = size(d,1)-1
            ni = 1
            ds = $TDataSeries{T,1}(nt, ni)
            copy!(ds.d, d)
            return ds
        end

        function $TDataSeries(d::AbstractArray{T,2}) where {T}
            nt = size(d,1)-1
            ni = size(d,2)
            ds = $TDataSeries{T,2}(nt, ni)
            copy!(ds.d, d)
            return ds
        end
    end
end

Base.:(==)(ds1::DataSeries, ds2::DataSeries) = (ds1.d == ds2.d)

function Base.show(io::IO, ds::DataSeries{T,N}) where {T,N}
    print(io, "DataSeries with data type ", T, " and ", N, " dimensions:\n")
    print(io, ds.d)
end


Base.parent(ds::DataSeries) = ds.d
Base.eltype(ds::DataSeries{T,N}) where {T,N} = T
Base.ndims(ds::DataSeries{T,N}) where {T,N} = N

GeometricBase.ntime(ds::DataSeries) = lastindex(ds.d, 1) - 1
GeometricBase.nsamples(ds::DataSeries) = size(ds.d, 2)

# errmsg_size(ds::DataSeries) = error("size not supported for data series with axes $(axes(ds))")
# Base.size(ds::DataSeries) = errmsg_size(ds)
# Base.size(ds::DataSeries, d) = errmsg_size(ds)
Base.size(ds::DataSeries) = size(ds.d)
Base.size(ds::DataSeries, d) = size(ds.d, d)


function Base.similar(ds::DS) where {T, DS <: DataSeries{T,1}}
    newds = DS(ntime(ds), nsamples(ds))
    for j in eachindex(ds)
        newds[j] = zero(ds[j])
    end
    return newds
end

function Base.similar(ds::DS) where {T, DS <: DataSeries{T,2}}
    newds = DS(ntime(ds), nsamples(ds))
    for k in axes(ds,2)
        for j in axes(ds,1)
            newds[j,k] = zero(ds[j,k])
        end
    end
    return newds
end

Base.eachindex(::IndexCartesian, ds::DataSeries) = CartesianIndices(axes(ds))
# Base.eachindex(::IndexLinear, ds::DataSeries) = axes(ds, 1)

Base.firstindex(ds::DataSeries{T,1}) where {T} = 0
Base.firstindex(ds::DataSeries{T,N}) where {T,N} = firstindex(ds.d)

Base.firstindex(ds::DataSeries{T,1}, d) where {T} = d ≥ 1 && d ≤ 1 ? 0 : 1
Base.firstindex(ds::DataSeries{T,2}, d) where {T} = d ≥ 1 && d ≤ 2 ? (0, 1)[d] : 1

Base.lastindex(ds::DataSeries{T,1}) where {T} = ntime(ds)
Base.lastindex(ds::DataSeries{T,N}) where {T,N} = lastindex(ds.d)

Base.lastindex(ds::DataSeries{T,1}, d) where {T} = d ≥ 1 && d ≤ 1 ? ntime(ds) : 1
Base.lastindex(ds::DataSeries{T,2}, d) where {T} = d ≥ 1 && d ≤ 2 ? (ntime(ds), nsamples(ds))[d] : 1

@inline Base.axes(ds::DataSeries{T,1}) where {T} = (0:ntime(ds),)
@inline Base.axes(ds::DataSeries{T,2}) where {T} = (0:ntime(ds), 1:nsamples(ds))
@inline Base.axes(ds::DataSeries{T,N}, d) where {T,N} = d ≥ 1 && d ≤ N ? axes(ds)[d] : (1:1)

Base.strides(ds::DataSeries) = strides(ds.d)


function get_data!(ds::DataSeries{T,1}, x::T, n::Int, k::Int=1) where {T}
    @assert n ≥ 0 && n ≤ ntime(ds)
    @assert k == 1
    x .= ds.d[n+1]
end

function get_data!(ds::DataSeries{T,2}, x::T, n::Int, k::Int) where {T}
    @assert n ≥ 0 && n ≤ ntime(ds)
    @assert k ≥ 1 && k ≤ nsamples(ds)
    x .= ds.d[n+1,k]
end

function get_data!(ds::DataSeries{T,2}, x::Vector{T}, n::Int) where {T}
    @assert length(x) == nsamples(ds)
    @assert n ≥ 0 && n ≤ ntime(ds)
    @inbounds for k in eachindex(x)
        x[k] = ds.d[n+1,k]
    end
end

function set_data!(ds::DataSeries{T,1}, x::T, n::Int, k::Int=1) where {T}
    @assert n ≥ 0 && n ≤ ntime(ds)
    @assert k == 1
    if T <: Number
        ds.d[n+1] = x
    elseif T <: AbstractArray
        ds.d[n+1] .= x
    end
end

function set_data!(ds::DataSeries{T,1}, x::Vector{T}, n::Int, k::Int=1) where {T}
    @assert nsamples(ds) == length(x) == 1
    set_data!(ds, x[begin], n, k)
end

function set_data!(ds::DataSeries{T,2}, x::T, n::Int, k::Int) where {T}
    @assert n ≥ 0 && n ≤ ntime(ds)
    if T <: Number
        ds.d[n+1,k] = x
    elseif T <: AbstractArray
        ds.d[n+1,k] .= x
    end
end

function set_data!(ds::DataSeries{T,2}, x::Vector{T}, n::Int) where {T}
    @assert length(x) == nsamples(ds)
    @assert n ≥ 0 && n ≤ ntime(ds)
    @inbounds for k in axes(ds.d, 2)
        if T isa Number
            ds.d[n+1,k] = x[k]
        elseif T isa AbstractArray
            ds.d[n+1,k] .= x[k]
        end
    end
end

function GeometricBase.reset!(ds::DataSeries{T,1}) where {T}
    @inbounds ds[begin] = ds[end]
end

function GeometricBase.reset!(ds::DataSeries{T,2}) where {T}
    @inbounds for k in axes(ds,2)
        ds[begin,k] = ds[end,k]
    end
end


@inline function Base.getindex(ds::DataSeries{T,1}, j::Union{UnitRange,Int}) where {T}
    getindex(ds.d, j.+1)
end

@inline function Base.getindex(ds::DataSeries{T,1}, i::Union{Int,CartesianIndex}, j::Union{Int,CartesianIndex}) where {T}
    ds[j][i]
end

@inline function Base.getindex(ds::DataSeries{T,1}, i::Union{Int,CartesianIndex}, J::AbstractRange{Int}) where {T}
    [ds[j][i] for j in J]
end

@inline function Base.getindex(ds::DataSeries{T,1}, i::Union{Int,CartesianIndex}, J::Colon) where {T}
    OffsetArray([ds[j][i] for j in axes(ds,1)], axes(ds,1))
end

@inline Base.getindex(ds::DataSeries{T,1}, I, J::Colon) where {T} = getindex(ds, I, axes(ds,1))


@inline function Base.getindex(ds::DataSeries{T,2}, j::Union{UnitRange,Int}, k::Union{UnitRange,Int}) where {T}
    getindex(ds.d, j.+1, k)
end

@inline function Base.getindex(ds::DataSeries{T,2}, i::Union{Int,CartesianIndex}, j::Union{Int,CartesianIndex}, k::Union{Int,CartesianIndex}) where {T}
    ds[j,k][i]
end

@inline function Base.getindex(ds::DataSeries{T,2}, i::Union{Int,CartesianIndex}, J::AbstractRange{Int}, K::AbstractRange{Int}) where {T}
    [ds[j,k][i] for j in J, k in K]
end

@inline function Base.getindex(ds::DataSeries{T,2}, i::Union{Int,CartesianIndex}, J::Colon, K::AbstractRange{Int}) where {T}
    OffsetArray([ds[j,k][i] for j in axes(ds,1), k in K], axes(ds,1), K)
end

@inline function Base.getindex(ds::DataSeries{T,2}, i::Union{Int,CartesianIndex}, J::AbstractRange{Int}, k::Union{Int,CartesianIndex}) where {T}
    [ds[j,k][i] for j in J]
end

@inline function Base.getindex(ds::DataSeries{T,2}, i::Union{Int,CartesianIndex}, J::Colon, k::Union{Int,CartesianIndex}) where {T}
    OffsetArray([ds[j,k][i] for j in axes(ds,1)], axes(ds,1))
end

@inline function Base.getindex(ds::DataSeries{T,2}, i::Union{Int,CartesianIndex}, j::Union{Int,CartesianIndex}, K::AbstractRange{Int}) where {T}
    [ds[j,k][i] for k in K]
end

@inline Base.getindex(ds::DataSeries{T,2}, J::Colon, K::Colon) where {T} = getindex(ds, axes(ds,1), axes(ds,2))
@inline Base.getindex(ds::DataSeries{T,2}, J, K::Colon) where {T} = getindex(ds, J, axes(ds,2))
@inline Base.getindex(ds::DataSeries{T,2}, J::Colon, K) where {T} = getindex(ds, axes(ds,1), K)

@inline Base.getindex(ds::DataSeries{T,2}, I, J::Colon, K::Colon) where {T} = getindex(ds, I, axes(ds,1), axes(ds,2))
@inline Base.getindex(ds::DataSeries{T,2}, I, J::Colon, K) where {T} = getindex(ds, I, axes(ds,1), K)
@inline Base.getindex(ds::DataSeries{T,2}, I, J, K::Colon) where {T} = getindex(ds, I, J, axes(ds,2))


@inline function Base.setindex!(ds::DataSeries{T,1}, x, j::Union{UnitRange,Int}) where {T}
    setindex!(ds.d, x, j.+1)
end

@inline function Base.setindex!(ds::DataSeries{T,2}, x, j::Union{UnitRange,Int}, k::Union{UnitRange,Int}) where {T}
    setindex!(ds.d, x, j.+1, k)
end


# 
# @inline function Base.getindex(ds::DataSeries{T,2}, I::Union{Int,CartesianIndex}, J::Union{Int,CartesianIndex,AbstractRange{Int}}, K::Union{Int,CartesianIndex,AbstractRange{Int}}) where {T}
#     x = zeros(eltype(ds.d[begin]), (1 for i in 1:length(I))..., length(J))
#     for j in J
#         x[Tuple(I)..., j] = ds[j][I]
#     end
# end
# 
# similar(storagetype, shape)
# 


function fromarray(::Type{dsType}, d::AbstractArray{T,1}, ni::Int=1) where {dsType <: DataSeries, T <: Number}
    @assert ni == 1
    return dsType(d)
end

function fromarray(::Type{dsType}, d::AbstractArray{T,2}, ni::Int=1) where {dsType <: DataSeries, T <: Number}
    @assert ni == 1 || ni == size(d,2)
    
    if ni == 1
        ds = [Vector{T}(d[:,n]) for n in axes(d,2)]
    else
        ds = Array(d)
    end

    return dsType(ds)
end

function fromarray(::Type{dsType}, d::AbstractArray{T,N}, ni::Int=1) where {dsType <: DataSeries, T <: Number, N}
    @assert ni == 1 || ni == size(d)[end]
    @assert N ≥ 3
    if ni == 1
        AT = Array{T,ndims(d)-1}
        ds = [AT(d[axes(d)[1:end-1]..., n]) for n in axes(d)[end]]
    else
        AT = Array{T,ndims(d)-2}
        ds = [AT(d[axes(d)[1:end-2]..., n, k]) for n in axes(d)[end-1], k in axes(d)[end]]
    end

    return dsType(ds)
end
