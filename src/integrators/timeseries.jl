
immutable Timeseries{T} <: AbstractArray{T,1}
    n::Int
    t::Vector{T}

    Δt::T
    step::Int

    function Timeseries(n, Δt, step)
        @assert T <: Real
        @assert n > 0
        @assert step > 0
        @assert n ≥ step

        t = zeros(T, n+1)
        new(n, t, Δt, step)
    end
end

function Timeseries{T}(n::Int, Δt::T, step::Int=1)
    return Timeseries{T}(n, Δt, step)
end

Base.eltype{T}(t::Timeseries{T}) = T
Base.ndims(t::Timeseries) = 1
Base.size(t::Timeseries) = size(t.t)
Base.length(t::Timeseries) = length(t.t)
Base.endof(t::Timeseries) = length(t)
Base.indices(t::Timeseries) = (0:t.n,)
Base.indices(t::Timeseries, d) = indices(t)[d]
Base.strides(t::Timeseries) = (1,)
Base.stride(t::Timeseries, d) = strides(t)[d]

@inline function Base.setindex!(ts::Timeseries, t, i::Int)
    @boundscheck checkbounds(ts.t, i+1)
    @inbounds setindex!(ts.t, t, i+1)
end

@inline function Base.getindex(ts::Timeseries, i::Int)
    @boundscheck checkbounds(ts.t, i+1)
    @inbounds t = getindex(ts.t, i+1)
    return t
end

function compute_timeseries!(ts::Timeseries)
    for n in 1:ts.n
        ts[n] = ts[0] + n * ts.step * ts.Δt
    end
end
