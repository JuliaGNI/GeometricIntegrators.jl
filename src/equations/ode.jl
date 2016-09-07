
# TODO Add momentum maps, energy, ...

immutable ODE{T} <: Equation{T}
    d::UInt
    f::Function
    q₀::Array{T, 1}

    function ODE(d, f, q₀)
        @assert d == length(q₀)
        @assert T == eltype(q₀)

        new(d, f, q₀)
    end
end

function ODE(d, f, q₀)
    ODE{eltype(q₀)}(d, f, q₀)
end


immutable PODE{T} <: Equation{T}
    d::UInt
    f::Function
    g::Function
    q₀::Array{T, 1}
    p₀::Array{T, 1}

    function PODE(d, f, g, q₀, p₀)
        @assert d == length(q₀) == length(p₀)
        @assert T == eltype(q₀) == eltype(p₀)

        new(d, f, g, q₀, p₀)
    end
end


function PODE(d, f, g, q₀, p₀)
    PODE{eltype(q₀)}(d, f, g, q₀, p₀)
end
