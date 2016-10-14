

immutable DAE{T} <: Equation{T}
    m::Int
    n::Int
    f::Function
    u::Function
    ϕ::Function
    t₀::T
    q₀::Array{T,1}
    λ₀::Array{T,1}

    function DAE(m, n, f, u, ϕ, t₀, q₀, λ₀)
        @assert m == length(q₀)
        @assert n == length(λ₀)
        @assert m ≥ n
        @assert T == eltype(q₀)
        @assert T == eltype(λ₀)

        new(m, n, f, u, ϕ, t₀, q₀, λ₀)
    end
end

function DAE{T}(m::Integer, n::Integer, f::Function, u::Function, ϕ::Function, t₀::Real, q₀::Vector{T}, λ₀::Vector{T})
    DAE{T}(m, n, f, u, ϕ, t₀, q₀, λ₀)
end

function DAE{T}(m::Integer, n::Integer, f::Function, u::Function, ϕ::Function, q₀::Vector{T}, λ₀::Vector{T})
    DAE{T}(m, n, f, u, ϕ, 0, q₀, λ₀)
end

function DAE{T}(f::Function, u::Function, ϕ::Function, t₀::Real, q₀::Vector{T}, λ₀::Vector{T})
    DAE{T}(length(q₀), length(λ₀), f, u, ϕ, t₀, q₀, λ₀)
end

function DAE{T}(f::Function, u::Function, ϕ::Function, q₀::Vector{T}, λ₀::Vector{T})
    DAE{T}(length(q₀), length(λ₀), f, u, ϕ, 0, q₀, λ₀)
end
