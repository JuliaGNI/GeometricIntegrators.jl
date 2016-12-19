
immutable TimeSeries{T} <: AbstractArray{T,1}
    n::Int
    t::Vector{T}

    Δt::T
    step::Int

    function TimeSeries(n, Δt, step)
        @assert T <: Real
        @assert n > 0
        @assert step > 0
        @assert n ≥ step

        t = zeros(T, n+1)
        new(n, t, Δt, step)
    end
end

function TimeSeries{T}(n::Int, Δt::T, step::Int=1)
    return TimeSeries{T}(n, Δt, step)
end

Base.eltype{T}(ts::TimeSeries{T}) = T
Base.ndims(ts::TimeSeries) = 1
Base.size(ts::TimeSeries) = size(ts.t)
Base.length(ts::TimeSeries) = length(ts.t)
Base.endof(ts::TimeSeries) = length(ts)
Base.indices(ts::TimeSeries) = (0:ts.n,)
Base.indices(ts::TimeSeries, d) = indices(ts)[d]
Base.strides(ts::TimeSeries) = (1,)
Base.stride(ts::TimeSeries, d) = strides(ts)[d]

@inline function Base.setindex!(ts::TimeSeries, t, i::Int)
    @boundscheck checkbounds(ts.t, i+1)
    @inbounds setindex!(ts.t, t, i+1)
end

@inline function Base.getindex(ts::TimeSeries, i::Int)
    @boundscheck checkbounds(ts.t, i+1)
    @inbounds t = getindex(ts.t, i+1)
    return t
end

function compute_timeseries!{T}(ts::TimeSeries{T}, t₀::T)
    ts[0] = t₀
    for n in 1:ts.n
        ts[n] = ts[0] + n * ts.step * ts.Δt
    end
end
