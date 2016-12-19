
immutable DataSeries{T,N} <: AbstractArray{T,N}
    nd::Int
    nt::Int
    ni::Int
    d::Array{T,N}
    step::Int

    function DataSeries(nd, nt, ni, step)
        @assert T <: Number
        @assert nd > 0
        @assert nt > 0
        @assert ni > 0
        @assert nt ≥ step

        @assert N ∈ (2,3)

        if N == 2
            d = zeros(T, nd, nt+1)
        elseif N == 3
            d = zeros(T, nd, nt+1, ni)
        end

        new(nd, nt, ni, d, step)
    end
end

function DataSeries(T, nd::Int, nt::Int, ni::Int, step::Int=1)
    ni == 1 ? N = 2 : N = 3
    return DataSeries{T,N}(nd, nt, ni, step)
end

Base.eltype{T,N}(ds::DataSeries{T,N}) = T
Base.ndims{T,N}(ds::DataSeries{T,N}) = N
Base.size(ds::DataSeries) = size(ds.d)
Base.length(ds::DataSeries) = length(ds.d)
Base.endof(ds::DataSeries) = length(ds)
Base.indices{T}(ds::DataSeries{T,2}) = (1:ds.nd, 0:ds.nt)
Base.indices{T}(ds::DataSeries{T,3}) = (1:ds.nd, 0:ds.nt, 1:ds.ni)
Base.strides(ds::DataSeries) = (strides(ds.d))

function set_initial_conditions!(ds::DataSeries, equation::ODE)
    set_initial_conditions!(ds, equation.t₀, equation.q₀)
end

function set_initial_conditions!{T}(ds::DataSeries{T,2}, x::Array{T,1})
    @assert size(ds, 1) == size(x, 1)
    @inbounds for i in 1:size(ds, 1)
        ds[i,0] = x[i]
    end
end

function set_initial_conditions!{T}(ds::DataSeries{T,3}, x::Array{T,2})
    @assert size(ds, 1) == size(x, 1)
    @assert size(ds, 3) == size(x, 2)
    @inbounds for k in 1:size(ds, 3)
        for i in 1:size(ds, 1)
            ds[i,0,k] = x[i,k]
        end
    end
end

function copy_solution!{T}(x::Vector{T}, ds::DataSeries{T,2}, n, k=1)
    if mod(n, ds.step) == 0
        j = div(n, ds.step)+1
        @assert length(x) == size(ds.d, 1)
        @assert j ≤ size(ds.d, 2)
        @assert k == 1
        @inbounds for i in 1:size(ds.d, 1)
            ds.d[i,j] = x[i]
        end
    end
end

function copy_solution!{T}(x::Vector{T}, ds::DataSeries{T,3}, n, k)
    if mod(n, ds.step) == 0
        j = div(n, ds.step)+1
        @assert length(x) == size(ds.d, 1)
        @assert j ≤ size(ds.d, 2)
        @assert k ≤ size(ds.d, 3)
        @inbounds for i in 1:size(ds.d, 1)
            ds.d[i,j,k] = x[i]
        end
    end
end

function reset!{T}(ds::DataSeries{T,2})
    for i in 1:size(ds,1)
        ds[i,0] = ds[i,end]
    end
end

function reset!{T}(ds::DataSeries{T,3})
    for k in 1:size(ds,3)
        for i in 1:size(ds,1)
            ds[i,0,k] = ds[i,end,k]
        end
    end
end

@inline function Base.getindex{T}(ds::DataSeries{T,2}, i::Int, j::Int)
    @boundscheck checkbounds(ds.d, i, j+1)
    @inbounds r = getindex(ds.d, i, j+1)
    return r
end

@inline function Base.getindex{T}(ds::DataSeries{T,2}, j::Int)
    @boundscheck checkbounds(ds.d, :, j+1)
    @inbounds r = getindex(ds.d, :, j+1)
    return r
end

@inline function Base.getindex{T}(ds::DataSeries{T,3}, i::Int, j::Int, k::Int)
    @boundscheck checkbounds(ds.d, i, j+1, k)
    @inbounds r = getindex(ds.d, i, j+1, k)
    return r
end

@inline function Base.getindex{T}(ds::DataSeries{T,3}, j::Int, k::Int)
    @boundscheck checkbounds(ds.d, :, j+1, k)
    @inbounds r = getindex(ds.d, :, j+1, k)
    return r
end

@inline function Base.getindex{T}(ds::DataSeries{T,3}, k::Int)
    @boundscheck checkbounds(ds.d, :, :, k)
    @inbounds r = getindex(ds.d, :, :, k)
    return r
end

@inline function Base.setindex!{T}(ds::DataSeries{T,2}, x, i::Int, j::Int)
    @assert length(x) == 1
    @boundscheck checkbounds(ds.d, i, j+1)
    @inbounds setindex!(ds.d, x, i, j+1)
end

@inline function Base.setindex!{T}(ds::DataSeries{T,2}, x, j::Int)
    @assert ndims(x) == 1
    @assert length(x) == size(ds.d, 1)
    @boundscheck checkbounds(ds.d, :, j+1)
    @inbounds setindex!(ds.d, x, :, j+1)
end

@inline function Base.setindex!{T}(ds::DataSeries{T,3}, x, i::Int, j::Int, k::Int)
    @assert length(x) == 1
    @boundscheck checkbounds(ds.d, i, j+1, k)
    @inbounds setindex!(ds.d, x, i, j+1, k)
end

@inline function Base.setindex!{T}(ds::DataSeries{T,3}, x, j::Int, k::Int)
    @assert ndims(x) == 1
    @assert length(x) == size(ds.d, 1)
    @boundscheck checkbounds(ds.d, :, j+1, k)
    @inbounds setindex!(ds.d, x, :, j+1, k)
end

@inline function Base.setindex!{T}(ds::DataSeries{T,3}, x, k::Int)
    @assert ndims(x) == 2
    @assert size(x,1) == size(ds.d, 1)
    @assert size(x,2) == size(ds.d, 2)
    @boundscheck checkbounds(ds.d, :, :, k)
    @inbounds setindex!(ds.d, x, :, :, k)
end
