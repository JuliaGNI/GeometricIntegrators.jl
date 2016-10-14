
"Special Partitioned Ordinary Differential Equation"
immutable SPODE{T} <: Equation{T}
    d::Int
    f::Function
    g::Function
    t₀::T
    q₀::Array{T, 1}
    p₀::Array{T, 1}

    function SPODE(d, f, g, t₀, q₀, p₀)
        @assert d == length(q₀) == length(p₀)
        @assert T == eltype(q₀) == eltype(p₀)

        new(d, f, g, t₀, q₀, p₀)
    end
end


function SPODE{T}(d::Integer, f::Function, g::Function, t₀::Real, q₀::Vector{T}, p₀::Vector{T})
    SPODE{T}(d, f, g, t₀, q₀, p₀)
end

function SPODE{T}(d::Integer, f::Function, g::Function, q₀::Vector{T}, p₀::Vector{T})
    SPODE{T}(d, f, g, 0, q₀, p₀)
end

function SPODE{T}(f::Function, g::Function, t₀::Real, q₀::Vector{T}, p₀::Vector{T})
    @assert length(q₀) == length(p₀)
    SPODE{T}(length(q₀), f, g, t₀, q₀, p₀)
end

function SPODE{T}(f::Function, g::Function, q₀::Vector{T}, p₀::Vector{T})
    @assert length(q₀) == length(p₀)
    SPODE{T}(length(q₀), f, g, 0, q₀, p₀)
end
