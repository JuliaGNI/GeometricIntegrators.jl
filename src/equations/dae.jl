
# TODO Add momentum maps, energy, ...

immutable DAE{T} <: Equation{T}
    m::UInt
    n::UInt
    f::Function
    u::Function
    ϕ::Function
    x0::Array{T,1}
    λ0::Array{T,1}

    function DAE(m, n, f, u, ϕ, x0, λ0)
        @assert m == length(x0) == length(f(x0)) == length(u(x0, λ0))
        @assert n == length(λ0) == length(ϕ(x0))
        @assert m ≥ n
        @assert eltype(x0) == eltype(f(x0)) == eltype(u(x0, λ0))
        @assert eltype(λ0) == eltype(ϕ(x0))

        new(m, n, f, u, ϕ, x0, λ0)
    end
end


immutable PDAE{T} <: Equation{T}
    m::UInt
    n::UInt
    f::Function
    g::Function
    u::Function
    v::Function
    ϕ::Function
    q0::Array{T,1}
    p0::Array{T,1}
    λ0::Array{T,1}

    function PDAE(m, n, f, g, u, v, ϕ, q0, p0, λ0)
        @assert m == length(q0) == length(f(q0, p0)) == length(u(q0, p0, λ0))
        @assert m == length(p0) == length(g(q0, p0)) == length(v(q0, p0, λ0))
        @assert n == length(λ0) == length(ϕ(q0, p0))
        @assert 2m ≥ n
        @assert eltype(q0) == eltype(f(q0, p0)) == eltype(u(q0, p0, λ0))
        @assert eltype(p0) == eltype(g(q0, p0)) == eltype(v(q0, p0, λ0))
        @assert eltype(λ0) == eltype(ϕ(q0, p0))

        new(m, n, f, g, u, v, ϕ, q0, p0, λ0)
    end
end
