

immutable ODE{T} <: Equation{T}
    d::Int
    f::Function
    t₀::T
    q₀::Array{T, 1}

    function ODE(d, f, t₀, q₀)
        @assert d == length(q₀)
        @assert T == eltype(q₀)

        new(d, f, t₀, q₀)
    end
end

function ODE{T}(f::Function, t₀::Real, q₀::Vector{T})
    ODE{T}(length(q₀), f, t₀, q₀)
end

function ODE{T}(f::Function, q₀::Vector{T})
    ODE(f, 0, q₀)
end
