
abstract type StochasticDataSeries{T,N} <: AbstractArray{T,N} end


#Type for holding a stochastic data series:
# nd   - dimension of the space of data
# nt+1 - number of time steps in the DataSeries including the initial condition
# ns   - number of sample paths of the stochastic process
# ni   - number of the initial conditions
# d    - the array holding the data
# T    - type of the held data
# N    - the number of dimensions of the array d
struct SStochasticDataSeries{T,N} <: StochasticDataSeries{T,N}
    nd::Int
    nt::Int
    ns::Int
    ni::Int
    d::Array{T,N}

    function SStochasticDataSeries{T,N}(nd, nt, ns, ni) where {T,N}
        @assert T <: Number
        @assert nd > 0
        @assert nt ≥ 0
        @assert ns ≥ 0
        @assert ni ≥ 0

        @assert N ∈ (1,2,3,4)

        if N == 1
            d = zeros(T, nt+1)
        elseif N == 2
            d = zeros(T, nd, nt+1)
        elseif N == 3
            ni==1 ? d = zeros(T, nd, nt+1, ns) : d = zeros(T, nd, nt+1, ni)
        elseif N==4
            d = zeros(T, nd, nt+1, ns, ni)
        end

        new(nd, nt, ns, ni, d)
    end

    function SStochasticDataSeries{T,N}(nd, nt, ns, ni, d) where {T,N}
        new(nd, nt, ns, ni, d)
    end
end



function SStochasticDataSeries(T, nd::Int, nt::Int, ns::Int, ni::Int)

    if nd==ns==ni==1
        N=1
    elseif ns==ni==1
        N=2
    elseif ns!=1!=ni
        N=4
    else
        N=3
    end

    return SStochasticDataSeries{T,N}(nd, nt, ns, ni)
end



function SStochasticDataSeries(d::Array{T,1}) where {T}
    nd = 1
    nt = size(d,1)-1
    ns = 1
    ni = 1
    return SStochasticDataSeries{T,1}(nd, nt, ns, ni, d)
end



function SStochasticDataSeries(d::Array{T,2}) where {T}
    nd = size(d,1)
    nt = size(d,2)-1
    ns = 1
    ni = 1
    return SStochasticDataSeries{T,2}(nd, nt, ns, ni, d)
end


# A 3-dimensional matrix d can represent a data series with ns=1 or ni=1.
# The second argument specifies whether there are multiple initial conditions
# (IC=true, so ni>1) or a single one (default IC=false, so ni=1)
function SStochasticDataSeries(d::Array{T,3}; IC=false) where {T}
    nd = size(d,1)
    nt = size(d,2)-1

    if IC==false
        ns = size(d,3)
        ni = 1
    else
        ni = size(d,3)
        ns = 1
    end

    return SStochasticDataSeries{T,3}(nd, nt, ns, ni, d)
end



function SStochasticDataSeries(d::Array{T,4}) where {T}
    nd = size(d,1)
    nt = size(d,2)-1
    ns = size(d,3)
    ni = size(d,4)
    return SStochasticDataSeries{T,4}(nd, nt, ns, ni, d)
end


#TODO: PStochasticDataSeries


function similar(ds::StochasticDataSeries{T,N}) where {T,N}
    typeof(ds){T,N}(ds.nd, ds.nt, ds.ns, ds.ni)
end

Base.eltype(ds::StochasticDataSeries{T,N}) where {T,N} = T
Base.ndims(ds::StochasticDataSeries{T,N}) where {T,N} = N
Base.size(ds::StochasticDataSeries) = size(ds.d)
Base.length(ds::StochasticDataSeries) = length(ds.d)
Base.endof(ds::StochasticDataSeries{T,1}) where {T} = (ds.nt)
Base.endof(ds::StochasticDataSeries{T,2}) where {T} = (ds.nd, ds.nt)
Base.endof(ds::StochasticDataSeries{T,3}) where {T} = ds.ni==1 ? (ds.nd, ds.nt, ds.ns) : (ds.nd, ds.nt, ds.ni)
Base.endof(ds::StochasticDataSeries{T,4}) where {T} = (ds.nd, ds.nt, ds.ns, ds.ni)
Base.indices(ds::StochasticDataSeries{T,1}) where {T} = (0:ds.nt)
Base.indices(ds::StochasticDataSeries{T,2}) where {T} = (1:ds.nd, 0:ds.nt)
Base.indices(ds::StochasticDataSeries{T,3}) where {T} = ds.ni==1 ? (1:ds.nd, 0:ds.nt, 1:ds.ns) : (1:ds.nd, 0:ds.nt, 1:ds.ni)
Base.indices(ds::StochasticDataSeries{T,4}) where {T} = (1:ds.nd, 0:ds.nt, 1:ds.ns, 1:ds.ni)
Base.strides(ds::StochasticDataSeries) = (strides(ds.d))





function get_data!(ds::StochasticDataSeries{T,1}, n) where {T}
    j = n+1
    @assert j ≤ size(ds.d, 1)
    @inbounds return ds.d[j]
end

function get_data!(ds::StochasticDataSeries{T,2}, x::Union{Array{T,1}, Array{Double{T},1}}, n) where {T}
    j = n+1
    @assert j ≤ size(ds.d, 2)
    @inbounds for i in eachindex(x)
        x[i] = ds.d[i,j]
    end
end

function get_data!(ds::StochasticDataSeries{T,3}, x::Union{Array{T,2}, Array{Double{T},2}}, n) where {T}
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

function get_data!(ds::StochasticDataSeries{T,3}, x::Union{Array{T,1}, Array{Double{T},1}}, n, k) where {T}
    j = n+1
    @assert size(x,1) == size(ds.d, 1)
    @assert j ≤ size(ds.d, 2)
    @assert k ≤ size(ds.d, 3)
    @inbounds for i in eachindex(x)
        x[i] = ds.d[i,j,k]
    end
end


function get_data!(ds::StochasticDataSeries{T,4}, x::Union{Array{T,3}, Array{Double{T},3}}, n) where {T}
    j = n+1
    @assert size(x,1) == size(ds.d, 1)
    @assert size(x,2) == size(ds.d, 3)
    @assert size(x,3) == size(ds.d, 4)
    @assert j ≤ size(ds.d, 2)
    @inbounds for m in 1:size(x,3)
        for k in 1:size(x, 2)
            for i in 1:size(x, 1)
                x[i,k,m] = ds.d[i,j,k,m]
            end
        end
    end
end


# Copies the data for the time step n and initial condition m from ds.d to x
function get_data!(ds::StochasticDataSeries{T,4}, x::Union{Array{T,2}, Array{Double{T},2}}, n, m) where {T}
    j = n+1
    @assert size(x,1) == size(ds.d, 1)
    @assert size(x,2) == size(ds.d, 3)
    @assert j ≤ size(ds.d, 2)
    @assert m ≤ size(ds.d, 4)
    @inbounds for k in 1:size(x, 2)
        for i in 1:size(x, 1)
            x[i,k] = ds.d[i,j,k,m]
        end
    end
end


function get_data!(ds::StochasticDataSeries{T,4}, x::Union{Array{T,1}, Array{Double{T},1}}, n, k, m) where {T}
    j = n+1
    @assert size(x,1) == size(ds.d, 1)
    @assert j ≤ size(ds.d, 2)
    @assert k ≤ size(ds.d, 3)
    @assert m ≤ size(ds.d, 4)
    @inbounds for i in eachindex(x)
        x[i] = ds.d[i,j,k,m]
    end
end



# Assigns a number to the nth time slot in a one-dimensional data series
function set_data!(ds::StochasticDataSeries{T,1}, x::Union{T, Double{T}}, n) where {T}
    j = n+1
    @assert j ≤ size(ds.d, 1)
    @inbounds ds.d[j] = x
end

# Assigns a number to the nth time slot in a one-dimensional data series,
# when the number is passed as a single-entry array
function set_data!(ds::StochasticDataSeries{T,1}, x::Union{Array{T,1}, Array{Double{T},1}}, n) where {T}
    @assert length(x)==1
    set_data!(ds,x[1],n)
end

# Assigns a vector to the nth time slot in a two-dimensional data series
function set_data!(ds::StochasticDataSeries{T,2}, x::Union{Array{T,1}, Array{Double{T},1}}, n) where {T}
    j = n+1
    @assert length(x) == size(ds.d, 1)
    @assert j ≤ size(ds.d, 2)
    @inbounds for i in 1:size(ds.d, 1)
        ds.d[i,j] = x[i]
    end
end

# Assigns a two dimensional array (nd x ns or nd x ni) to the nth time slot in a 3-dimensional data series
function set_data!(ds::StochasticDataSeries{T,3}, x::Union{Array{T,2}, Array{Double{T},2}}, n) where {T}
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

# Assigns a one dimensional array (nd elements) to the nth time slot and the kth sample path in a 3-dimensional data series
function set_data!(ds::StochasticDataSeries{T,3}, x::Union{Array{T,1}, Array{Double{T},1}}, n, k) where {T}
    j = n+1
    @assert size(x,1) == size(ds.d, 1)
    @assert j ≤ size(ds.d, 2)
    @assert k ≤ size(ds.d, 3)
    @inbounds for i in 1:size(ds.d, 1)
        ds.d[i,j,k] = x[i]
    end
end

# Assigns a one dimensional array (nd elements) to the nth time slot and to all sample paths in a 3-dimensional data series
function set_data!(ds::StochasticDataSeries{T,3}, x::Union{Array{T,1}, Array{Double{T},1}}, n) where {T}
    j = n+1
    @assert size(x,1) == size(ds.d, 1)
    @assert j ≤ size(ds.d, 2)
    @inbounds for k in 1:size(ds.d, 3)
        @inbounds for i in 1:size(ds.d, 1)
            ds.d[i,j,k] = x[i]
        end
    end
end

# Assign an nd x ns x ni array to the nth time slot in a 4-dimensional data series
function set_data!(ds::StochasticDataSeries{T,4}, x::Union{Array{T,3}, Array{Double{T},3}}, n) where {T}
    j = n+1
    @assert size(x,1) == size(ds.d, 1)
    @assert size(x,2) == size(ds.d, 3)
    @assert size(x,3) == size(ds.d, 4)
    @assert j ≤ size(ds.d, 2)
    @inbounds for m in 1:size(x,3)
        for k in 1:size(x, 2)
            for i in 1:size(x, 1)
                ds.d[i,j,k,m] = x[i,k,m]
            end
        end
    end
end


# Assings the data stored in x to the data series ds as the time step n and for
# the initial condition m
function set_data!(ds::StochasticDataSeries{T,4}, x::Union{Array{T,2}, Array{Double{T},2}}, n, m) where {T}
    j = n+1
    @assert size(x,1) == size(ds.d, 1)
    @assert size(x,2) == size(ds.d, 3)
    @assert j ≤ size(ds.d, 2)
    @assert m ≤ size(ds.d, 4)
    @inbounds for k in 1:size(x, 2)
        for i in 1:size(x, 1)
            ds.d[i,j,k,m] = x[i,k]
        end
    end
end


# Assings an nd x ni matrix to the data series ds as the time step n and for
# all sample paths
function set_data!(ds::StochasticDataSeries{T,4}, x::Union{Array{T,2}, Array{Double{T},2}}, n) where {T}
    j = n+1
    @assert size(x,1) == size(ds.d, 1)
    @assert size(x,2) == size(ds.d, 4)
    @assert j ≤ size(ds.d, 2)
    @inbounds for m in 1:size(ds.d, 3)
        @inbounds for k in 1:size(x, 2)
            for i in 1:size(x, 1)
                ds.d[i,j,m,k] = x[i,k]
            end
        end
    end
end


function set_data!(ds::StochasticDataSeries{T,4}, x::Union{Array{T,1}, Array{Double{T},1}}, n, k, m) where {T}
    j = n+1
    @assert size(x,1) == size(ds.d, 1)
    @assert j ≤ size(ds.d, 2)
    @assert k ≤ size(ds.d, 3)
    @assert m ≤ size(ds.d, 4)
    @inbounds for i in eachindex(x)
        ds.d[i,j,k,m] = x[i]
    end
end


function reset!(ds::StochasticDataSeries{T,1}) where {T}
    # ds[0]=ds[end-1] was wrong for some reason
    # ds[0]=ds[end] would be ok, but to be safe below using ds.d
    @inbounds ds.d[1] = ds.d[end]
end

function reset!(ds::StochasticDataSeries{T,2}) where {T}
    # when 'end' is passed to getindex, because of the shift it becomes end+1
    # in ds.d, therefore end-1 has to be passed to getindex instead
    @inbounds for i in 1:size(ds,1)
        ds[i,0] = ds[i,end-1]
    end
end

function reset!(ds::StochasticDataSeries{T,3}) where {T}
    @inbounds for k in 1:size(ds,3)
        for i in 1:size(ds,1)
            ds[i,0,k] = ds[i,end-1,k]
        end
    end
end

function reset!(ds::StochasticDataSeries{T,4}) where {T}
    @inbounds for m in 1:size(ds,4)
        for k in 1:size(ds,3)
            for i in 1:size(ds,1)
                ds[i,0,k,m] = ds[i,end-1,k,m]
            end
        end
    end
end

@inline function Base.getindex(ds::StochasticDataSeries{T,1}, j::Int) where {T}
    @boundscheck checkbounds(ds.d, j+1)
    @inbounds r = getindex(ds.d, j+1)
    return r
end

@inline function Base.getindex(ds::StochasticDataSeries{T,2}, i::Int, j::Int) where {T}
    @boundscheck checkbounds(ds.d, i, j+1)
    @inbounds r = getindex(ds.d, i, j+1)
    return r
end

@inline function Base.getindex(ds::StochasticDataSeries{T,2}, j::Int) where {T}
    @boundscheck checkbounds(ds.d, :, j+1)
    @inbounds r = getindex(ds.d, :, j+1)
    return r
end

@inline function Base.getindex(ds::StochasticDataSeries{T,3}, i::Int, j::Int, k::Int) where {T}
    @boundscheck checkbounds(ds.d, i, j+1, k)
    @inbounds r = getindex(ds.d, i, j+1, k)
    return r
end

@inline function Base.getindex(ds::StochasticDataSeries{T,3}, j::Int, k::Int) where {T}
    @boundscheck checkbounds(ds.d, :, j+1, k)
    @inbounds r = getindex(ds.d, :, j+1, k)
    return r
end

@inline function Base.getindex(ds::StochasticDataSeries{T,3}, k::Int) where {T}
    @boundscheck checkbounds(ds.d, :, :, k)
    @inbounds r = getindex(ds.d, :, :, k)
    return r
end

@inline function Base.getindex(ds::StochasticDataSeries{T,4}, i::Int, j::Int, k::Int, m::Int) where {T}
    @boundscheck checkbounds(ds.d, i, j+1, k, m)
    @inbounds r = getindex(ds.d, i, j+1, k, m)
    return r
end

@inline function Base.getindex(ds::StochasticDataSeries{T,4}, j::Int, k::Int, m::Int) where {T}
    @boundscheck checkbounds(ds.d, :, j+1, k, m)
    @inbounds r = getindex(ds.d, :, j+1, k, m)
    return r
end

@inline function Base.getindex(ds::StochasticDataSeries{T,4}, j::Int, m::Int) where {T}
    @boundscheck checkbounds(ds.d, :, j+1, :, m)
    @inbounds r = getindex(ds.d, :, j+1, :, m)
    return r
end

@inline function Base.getindex(ds::StochasticDataSeries{T,4}, m::Int) where {T}
    @boundscheck checkbounds(ds.d, :, :, :, m)
    @inbounds r = getindex(ds.d, :, :, :, m)
    return r
end

@inline function Base.setindex!(ds::StochasticDataSeries{T,1}, x, j::Int) where {T}
    @assert length(x) == 1
    @boundscheck checkbounds(ds.d, j+1)
    @inbounds setindex!(ds.d, x, j+1)
end

@inline function Base.setindex!(ds::StochasticDataSeries{T,2}, x, i::Int, j::Int) where {T}
    @assert length(x) == 1
    @boundscheck checkbounds(ds.d, i, j+1)
    @inbounds setindex!(ds.d, x, i, j+1)
end

@inline function Base.setindex!(ds::StochasticDataSeries{T,2}, x, j::Int) where {T}
    @assert ndims(x) == 1
    @assert length(x) == size(ds.d, 1)
    @boundscheck checkbounds(ds.d, :, j+1)
    @inbounds setindex!(ds.d, x, :, j+1)
end

@inline function Base.setindex!(ds::StochasticDataSeries{T,3}, x, i::Int, j::Int, k::Int) where {T}
    @assert length(x) == 1
    @boundscheck checkbounds(ds.d, i, j+1, k)
    @inbounds setindex!(ds.d, x, i, j+1, k)
end

@inline function Base.setindex!(ds::StochasticDataSeries{T,3}, x, j::Int, k::Int) where {T}
    @assert ndims(x) == 1
    @assert length(x) == size(ds.d, 1)
    @boundscheck checkbounds(ds.d, :, j+1, k)
    @inbounds setindex!(ds.d, x, :, j+1, k)
end

@inline function Base.setindex!(ds::StochasticDataSeries{T,3}, x, j::Int) where {T}
    @assert ndims(x) == 2
    @assert size(x,1) == size(ds.d, 1)
    @assert size(x,2) == size(ds.d, 3)
    @boundscheck checkbounds(ds.d, :, j+1, :)
    @inbounds setindex!(ds.d, x, :, j+1, :)
end


@inline function Base.setindex!(ds::StochasticDataSeries{T,4}, x, i::Int, j::Int, k::Int, m::Int) where {T}
    @assert length(x) == 1
    @boundscheck checkbounds(ds.d, i, j+1, k, m)
    @inbounds setindex!(ds.d, x, i, j+1, k, m)
end

@inline function Base.setindex!(ds::StochasticDataSeries{T,4}, x, j::Int, k::Int, m::Int) where {T}
    @assert ndims(x) == 1
    @assert length(x) == size(ds.d, 1)
    @boundscheck checkbounds(ds.d, :, j+1, k, m)
    @inbounds setindex!(ds.d, x, :, j+1, k, m)
end


@inline function Base.setindex!(ds::StochasticDataSeries{T,4}, x, j::Int, m::Int) where {T}
    @assert ndims(x) == 2
    @assert size(x,1) == size(ds.d, 1)
    @assert size(x,2) == size(ds.d, 3)
    @boundscheck checkbounds(ds.d, :, j+1, :, m)
    @inbounds setindex!(ds.d, x, :, j+1, :, m)
end


@inline function Base.setindex!(ds::StochasticDataSeries{T,4}, x, j::Int) where {T}
    @assert ndims(x) == 3
    @assert size(x,1) == size(ds.d, 1)
    @assert size(x,2) == size(ds.d, 3)
    @assert size(x,3) == size(ds.d, 4)
    @boundscheck checkbounds(ds.d, :, j+1, :, :)
    @inbounds setindex!(ds.d, x, :, j+1, :, :)
end

# TODO Implement convert() to/from array.
