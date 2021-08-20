
struct DataSeries{T <: AbstractData, N} <: AbstractDataSeries{T, N, Array{T,N}}
    d::Array{T,N}

    function DataSeries{T,N}(nt, ni) where {T <: AbstractData, N}
        @assert nt ≥ 0
        @assert ni > 0

        @assert N ∈ (1,2)

        if N == 1
            d = Array{T,N}(undef, nt+1)
        elseif N == 2
            d = Array{T,N}(undef, nt+1, ni)
        end

        new(d)
    end
end

Base.parent(ds::DataSeries) = ds.d

DataSeries(::Type{T}, nt::Int, ni::Int=1) where {T} = DataSeriesConstructor(DataSeries, T, nt, ni)
DataSeries(::Type{T}, nd::Int, nt::Int, ni::Int) where {T} = DataSeriesConstructor(DataSeries, T, nd, nt, ni)
DataSeries(q₀::T, nt::Int, ni::Int=1) where {T <: AbstractData} = DataSeriesConstructor(DataSeries, q₀, nt, ni)
DataSeries(q₀::Vector{T}, nt::Int, ni::Int=1) where {T <: Number} = DataSeriesConstructor(DataSeries, q₀, nt, ni)
DataSeries(q₀::Vector{T}, nt::Int, ni::Int=length(q₀)) where {T <: AbstractArray{<:Number}} = DataSeriesConstructor(DataSeries, q₀, nt, ni)
DataSeries(d::AbstractArray) = DataSeriesConstructor(DataSeries, d)
