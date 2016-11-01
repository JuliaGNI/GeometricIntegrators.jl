
"Partitioned Ordinary Differential Equation"
immutable PODE{T} <: Equation{T}
    d::Int
    f::Function
    g::Function
    t₀::T
    q₀::Array{T, 1}
    p₀::Array{T, 1}

    function PODE(d, f, g, t₀, q₀, p₀)
        @assert d == length(q₀) == length(p₀)
        @assert T == eltype(q₀) == eltype(p₀)

        new(d, f, g, t₀, q₀, p₀)
    end
end


function PODE{T}(f::Function, g::Function, t₀::Real, q₀::Vector{T}, p₀::Vector{T})
    @assert length(q₀) == length(p₀)
    PODE{T}(length(q₀), f, g, t₀, q₀, p₀)
end

function PODE{T}(f::Function, g::Function, q₀::Vector{T}, p₀::Vector{T})
    PODE(f, g, 0, q₀, p₀)
end
