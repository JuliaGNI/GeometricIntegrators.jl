
abstract Solution{T} <: DenseArray{T,2}

function Solution(equation::ODE, Δt::Real, ntime::Int, nsave::Int=1)
    SolutionODE(equation, Δt, ntime, nsave)
end

function Solution(equation::PODE, Δt::Real, ntime::Int, nsave::Int=1)
    SolutionPODE(equation, Δt, ntime, nsave)
end

function Solution(equation::DAE, Δt::Real, ntime::Int, nsave::Int=1)
    SolutionDAE(equation, Δt, ntime, nsave)
end

function Solution(equation::PDAE, Δt::Real, ntime::Int, nsave::Int=1)
    SolutionPDAE(equation, Δt, ntime, nsave)
end

function Solution(equation::Equation, Δt::Real, ntime::Int, nsave::Int=1)
    error("No solution found for equation ", equation)
end


Base.size(s::Solution) = size(s.x)

Base.indices(s::Solution) = (1:s.d, 0:s.n)

function Base.indices(s::Solution, d)
    if d == 1
        return 1:s.d
    elseif d == 2
        return 0:s.n
    end
end

function Base.stride(s::Solution, d)
    if d == 1
        return 1
    elseif d == 2
        return s.d
    end
end


@inline function Base.getindex(s::Solution, i::Int, j::Int)
    @boundscheck checkbounds(s.x, i, j+1)
    @inbounds r = getindex(s.x, i, j+1)
    return r
end

@inline function Base.setindex!(s::Solution, x, i::Int, j::Int)
    @boundscheck checkbounds(s.x, i, j+1)
    @inbounds setindex!(s.x, x, i, j+1)
end

# TODO Implement similar() and convert() to/from array functions.

# TODO Add solver status information to all Solutions.

immutable SolutionODE{T} <: Solution{T}
    d::Int
    n::Int
    x::Array{T,2}
    Δt::T
    ntime::Int
    nsave::Int

    function SolutionODE(d, Δt, ntime, nsave)
        @assert T <: Real
        @assert d > 0
        @assert nsave > 0
        @assert ntime ≥ nsave
        @assert mod(ntime, nsave) == 0

        n = div(ntime, nsave)
        new(d, n, zeros(T, d, n+1), Δt, ntime, nsave)
    end
end

function SolutionODE(equation::ODE, Δt::Real, ntime::Int, nsave::Int=1)
    T = eltype(equation.x0)
    s = SolutionODE{T}(equation.d, Δt, ntime, nsave)
    setInitialConditions(s, equation)
    return s
end

function setInitialConditions(solution::SolutionODE, equation::ODE)
    solution[1:equation.d, 0] = equation.x0
end

@inline function Base.getindex(s::SolutionODE, j::Int)
    @boundscheck checkbounds(s.x, 1:s.d, j+1)
    @inbounds r = getindex(s.x, 1:s.d, j+1)
    return r
end

@inline function Base.setindex!(s::SolutionODE, x, j::Int)
    @assert length(x) == s.d
    @boundscheck checkbounds(s.x, 1:s.d, j+1)
    @inbounds setindex!(s.x, x, 1:s.d, j+1)
end

function reset(s::SolutionODE)
    # TODO
end


immutable SolutionPODE{T} <: Solution{T}
    d::Int
    n::Int
    q::Array{T,2}
    p::Array{T,2}
    Δt::T
    ntime::Int
    nsave::Int

    function SolutionPODE(d, Δt, ntime, nsave)
        @assert T <: Real
        @assert d > 0
        @assert nsave > 0
        @assert ntime ≥ nsave
        @assert mod(ntime, nsave) == 0

        n = div(ntime, nsave)
        new(d, n, zeros(T, d, n+1), zeros(T, d, n+1), Δt, ntime, nsave)
    end
end

function SolutionPODE(equation::PODE, Δt::Real, ntime::Int, nsave::Int=1)
    T1 = eltype(equation.q0)
    T2 = eltype(equation.p0)
    @assert T1 == T2
    SolutionPODE{T1}(equation.d, Δt, ntime, nsave)
end

function Base.getindex(s::SolutionPODE, i::Int)
    # TODO
end

function reset(s::SolutionPODE)
    # TODO
end


immutable SolutionDAE{T} <: Solution{T}
    d::Int
    m::Int
    n::Int
    x::Array{T,2}
    λ::Array{T,2}
    Δt::T
    ntime::Int
    nsave::Int

    function SolutionDAE(d, m, Δt, ntime, nsave)
        @assert T <: Real
        @assert d > 0
        @assert m > 0
        @assert nsave > 0
        @assert ntime ≥ nsave
        @assert mod(ntime, nsave) == 0

        n = div(ntime, nsave)
        new(d, m, n, zeros(T, d, n), zeros(T, m, n), Δt, ntime, nsave)
    end
end

function SolutionDAE(equation::DAE, Δt::Real, ntime::Int, nsave::Int=1)
    T1 = eltype(equation.x0)
    T2 = eltype(equation.λ0)
    @assert T1 == T2
    SolutionDAE{T1}(equation.m, equation.n, Δt, ntime, nsave)
end

function Base.getindex(s::SolutionDAE, i::Int)
    # TODO
end

function reset(s::SolutionDAE)
    # TODO
end


immutable SolutionPDAE{T} <: Solution{T}
    d::Int
    m::Int
    n::Int
    q::Array{T,2}
    p::Array{T,2}
    λ::Array{T,2}
    Δt::T
    ntime::Int
    nsave::Int

    function SolutionPDAE(d, m, Δt, ntime, nsave)
        @assert T <: Real
        @assert d > 0
        @assert m > 0
        @assert nsave > 0
        @assert ntime ≥ nsave
        @assert mod(ntime, nsave) == 0

        n = div(ntime, nsave)
        new(d, m, n, zeros(T, d, n+1), zeros(T, d, n+1), zeros(T, m, n+1), Δt, ntime, nsave)
    end
end

function SolutionPDAE(equation::PDAE, Δt::Real, ntime::Int, nsave::Int=1)
    T1 = eltype(equation.q0)
    T2 = eltype(equation.p0)
    T3 = eltype(equation.λ0)
    @assert T1 == T2 == T3
    SolutionPDAE{T1}(equation.m, equation.n, Δt, ntime, nsave)
end

function Base.getindex(s::SolutionPDAE, i::Int)
    # TODO
end

function reset(s::SolutionPDAE)
    # TODO
end
