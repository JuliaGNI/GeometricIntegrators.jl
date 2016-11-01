
"Partitioned Differential Algebraic Equation"
immutable PDAE{T} <: Equation{T}
    m::Int
    n::Int
    f::Function
    g::Function
    u::Function
    v::Function
    ϕ::Function
    t₀::T
    q₀::Array{T,1}
    p₀::Array{T,1}
    λ₀::Array{T,1}

    function PDAE(m, n, f, g, u, v, ϕ, t₀, q₀, p₀, λ₀)
        @assert m == length(q₀)
        @assert m == length(p₀)
        @assert n == length(λ₀)
        @assert 2m ≥ n
        @assert T == eltype(q₀)
        @assert T == eltype(p₀)
        @assert T == eltype(λ₀)

        new(m, n, f, g, u, v, ϕ, t₀, q₀, p₀, λ₀)
    end
end

function PDAE{T}(f::Function, g::Function, u::Function, v::Function, ϕ::Function,
                 t₀::Real, q₀::Vector{T}, p₀::Vector{T}, λ₀::Vector{T})
    @assert length(q₀) == length(p₀)
    PDAE{T}(length(q₀), length(λ₀), f, g, u, v, ϕ, t₀, q₀, p₀, λ₀)
end

function PDAE{T}(f::Function, g::Function, u::Function, v::Function, ϕ::Function,
                 q₀::Vector{T}, p₀::Vector{T}, λ₀::Vector{T})
    PDAE(f, g, u, v, ϕ, 0, q₀, p₀, λ₀)
end
