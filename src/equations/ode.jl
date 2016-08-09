
# TODO Add momentum maps, energy, ...

immutable ODE <: Equation
    d::UInt
    f::Function
    x0::Array{Real, 1}

    function ODE(d, f, x0)
        @assert d == length(x0) == length(f(x0))
        @assert eltype(x0) == eltype(f(x0))

        new(d, f, x0)
    end
end


immutable PODE <: Equation
    d::UInt
    f::Function
    g::Function
    q0::Array{Real, 1}
    p0::Array{Real, 1}

    function PODE(d, f, g, q0, p0)
        @assert d == length(q0) == length(f(q0, p0))
        @assert d == length(p0) == length(g(q0, p0))
        @assert eltype(q0) == eltype(f(q0, p0))
        @assert eltype(p0) == eltype(g(q0, p0))

        new(d, f, g, q0, p0)
    end
end
