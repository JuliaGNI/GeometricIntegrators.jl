
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

Base.eltype{T}(ts::Timeseries{T}) = T
Base.ndims(ts::Timeseries) = 1
Base.size(ts::Timeseries) = size(ts.t)
Base.length(ts::Timeseries) = length(ts.t)
Base.endof(ts::Timeseries) = length(ts)
Base.indices(ts::Timeseries) = (0:ts.n,)
Base.indices(ts::Timeseries, d) = indices(ts)[d]
Base.strides(ts::Timeseries) = (1,)
Base.stride(ts::Timeseries, d) = strides(ts)[d]

@inline function Base.setindex!(ts::Timeseries, t, i::Int)
    @boundscheck checkbounds(ts.t, i+1)
    @inbounds setindex!(ts.t, t, i+1)
end

@inline function Base.getindex(ts::Timeseries, i::Int)
    @boundscheck checkbounds(ts.t, i+1)
    @inbounds t = getindex(ts.t, i+1)
    return t
end

function compute_timeseries!{T}(ts::Timeseries{T}, t₀::T)
    ts[0] = t₀
    for n in 1:ts.n
        ts[n] = ts[0] + n * ts.step * ts.Δt
    end
end
