
abstract Solution{T}

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

Base.length(s::Solution) = div(s.ntime, s.nsave)
Base.endof(s::Solution) = length(s)


# TODO Add solver status information to all Solutions.

type SolutionODE{T} <: Solution{T}
    d::UInt
    x::Array{T,2}
    Δt::Real
    ntime::UInt
    nsave::UInt

    function SolutionODE(d, Δt, ntime, nsave)
        @assert T <: Real
        @assert ntime ≥ nsave
        @assert mod(ntime, nsave) == 0

        new(d, zeros(T, d, div(ntime, nsave)), Δt, ntime, nsave)
    end
end

function SolutionODE(equation::ODE, Δt::Real, ntime::Int, nsave::Int=1)
    T = eltype(equation.x0)
    SolutionODE{T}(equation.d, Δt, ntime, nsave)
end

function Base.getindex(s::SolutionODE, i::Int)
    # TODO
end

function reset(s::SolutionODE)
    # TODO
end


type SolutionPODE{T} <: Solution{T}
    d::UInt
    q::Array{T,2}
    p::Array{T,2}
    Δt::Real
    ntime::UInt
    nsave::UInt

    function SolutionPODE(d, Δt, ntime, nsave)
        @assert T <: Real
        @assert ntime ≥ nsave
        @assert mod(ntime, nsave) == 0

        new(d, zeros(T, d, div(ntime, nsave)), zeros(T, d, div(ntime, nsave)), Δt, ntime, nsave)
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


type SolutionDAE{T} <: Solution{T}
    m::UInt
    n::UInt
    x::Array{T,2}
    λ::Array{T,2}
    Δt::Real
    ntime::UInt
    nsave::UInt

    function SolutionDAE(m, n, Δt, ntime, nsave)
        @assert T <: Real
        @assert ntime ≥ nsave
        @assert mod(ntime, nsave) == 0

        new(m, n, zeros(T, m, div(ntime, nsave)), zeros(T, n, div(ntime, nsave)), Δt, ntime, nsave)
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


type SolutionPDAE{T} <: Solution{T}
    m::UInt
    n::UInt
    q::Array{T,2}
    p::Array{T,2}
    λ::Array{T,2}
    Δt::Real
    ntime::UInt
    nsave::UInt

    function SolutionPDAE(m, n, Δt, ntime, nsave)
        @assert T <: Real
        @assert ntime ≥ nsave
        @assert mod(ntime, nsave) == 0

        new(m, n, zeros(T, m, div(ntime, nsave)), zeros(T, m, div(ntime, nsave)), zeros(T, n, div(ntime, nsave)), Δt, ntime, nsave)
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
