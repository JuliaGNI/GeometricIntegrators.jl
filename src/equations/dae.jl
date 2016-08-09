
# TODO Add momentum maps, energy, ...

immutable DAE <: Equation
    m::UInt
    n::UInt
    f::Function
    u::Function
    ϕ::Function
    x0::AbstractArray{Real, 1}
    λ0::AbstractArray{Real, 1}

    function DAE(m, n, f, u, ϕ, x0, λ0)
        @assert m == length(x0) == length(f(x0)) == length(u(x0, λ0))
        @assert n == length(λ0) == length(ϕ(x0))
        @assert m ≥ n

        new(m, n, f, u, ϕ, x0, λ0)
    end
end


immutable PDAE <: Equation
    m::UInt
    n::UInt
    f::Function
    g::Function
    u::Function
    v::Function
    ϕ::Function
    q0::AbstractArray{Real, 1}
    p0::AbstractArray{Real, 1}
    λ0::AbstractArray{Real, 1}

    function PDAE(m, n, f, g, u, v, ϕ, q0, p0, λ0)
        @assert m == length(q0) == length(f(q0, p0)) == length(u(q0, p0, λ0))
        @assert m == length(p0) == length(g(q0, p0)) == length(v(q0, p0, λ0))
        @assert n == length(λ0) == length(ϕ(q0, p0))
        @assert 2m ≥ n

        new(m, n, f, g, u, v, ϕ, q0, p0, λ0)
    end
end
