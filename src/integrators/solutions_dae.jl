
"Solution of a differential algebraic equation."
immutable SolutionDAE{T} <: Solution{T,3}
    d::Int
    m::Int
    n::Int
    t::Array{T,1}
    x::Array{T,2}
    λ::Array{T,2}
    ntime::Int
    nsave::Int

    function SolutionDAE(d, m, ntime, nsave)
        @assert T <: Real
        @assert d > 0
        @assert m > 0
        @assert nsave > 0
        @assert ntime ≥ nsave
        @assert mod(ntime, nsave) == 0

        n = div(ntime, nsave)
        t = zeros(T, n+1)
        x = zeros(T, d, n+1)
        λ = zeros(T, m, n+1)
        new(d, m, n, t, x, λ, ntime, nsave)
    end
end

function SolutionDAE(equation::DAE, ntime::Int, nsave::Int=1)
    T1 = eltype(equation.q₀)
    T2 = eltype(equation.λ₀)
    @assert T1 == T2
    SolutionDAE{T1}(equation.m, equation.n, ntime, nsave)
end

function Base.getindex(s::SolutionDAE, i::Int)
    # TODO
end

function reset(s::SolutionDAE)
    # TODO
end
