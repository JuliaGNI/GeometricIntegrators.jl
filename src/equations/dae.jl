
# TODO Add momentum maps, energy, ...

immutable DAE{T} <: Equation{T}
    m::Int
    n::Int
    f::Function
    u::Function
    ϕ::Function
    q₀::Array{T,1}
    λ₀::Array{T,1}

    function DAE(m, n, f, u, ϕ, q₀, λ₀)
        @assert m == length(q₀)
        @assert n == length(λ₀)
        @assert m ≥ n
        @assert T == eltype(q₀)
        @assert T == eltype(λ₀)

        new(m, n, f, u, ϕ, q₀, λ₀)
    end
end

function DAE(m, n, f, u, ϕ, q₀, λ₀)
    DAE{eltype(q₀)}(m, n, f, u, ϕ, q₀, λ₀)
end


immutable PDAE{T} <: Equation{T}
    m::Int
    n::Int
    f::Function
    g::Function
    u::Function
    v::Function
    ϕ::Function
    q₀::Array{T,1}
    p₀::Array{T,1}
    λ₀::Array{T,1}

    function PDAE(m, n, f, g, u, v, ϕ, q₀, p₀, λ₀)
        @assert m == length(q₀)
        @assert m == length(p₀)
        @assert n == length(λ₀)
        @assert 2m ≥ n
        @assert T == eltype(q₀)
        @assert T == eltype(p₀)
        @assert T == eltype(λ₀)

        new(m, n, f, g, u, v, ϕ, q₀, p₀, λ₀)
    end
end

function PDAE(m, n, f, g, u, v, ϕ, q₀, p₀, λ₀)
    PDAE{eltype(q₀)}(m, n, f, g, u, v, ϕ, q₀, p₀, λ₀)
end
