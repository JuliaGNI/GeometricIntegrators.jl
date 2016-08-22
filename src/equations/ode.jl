
# TODO Add momentum maps, energy, ...

immutable ODE{T} <: Equation{T}
    d::UInt
    f::Function
    x0::Array{T, 1}

    function ODE(d, f, x0)
        @assert d == length(x0)# == length(f(x0))
        @assert T == eltype(x0)# == eltype(f(x0))

        new(d, f, x0)
    end
end

function ODE(d, f, x0)
    ODE{eltype(x0)}(d, f, x0)
end


immutable PODE{T} <: Equation{T}
    d::UInt
    f::Function
    g::Function
    q0::Array{T, 1}
    p0::Array{T, 1}

    function PODE(d, f, g, q0, p0)
        @assert d == length(q0)# == length(f(q0, p0))
        @assert d == length(p0)# == length(g(q0, p0))
        @assert T == eltype(q0)# == eltype(f(q0, p0))
        @assert T == eltype(p0)# == eltype(g(q0, p0))

        new(d, f, g, q0, p0)
    end
end


function PODE(d, f, g, q0, p0)
    PODE{eltype(q0)}(d, f, g, q0, p0)
end
