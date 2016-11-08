
"Solution of a partitioned differential algebraic equation."
immutable SolutionPDAE{T} <: Solution{T,3}
    d::Int
    m::Int
    n::Int
    t::Array{T,1}
    q::Array{T,2}
    p::Array{T,2}
    λ::Array{T,2}
    ntime::Int
    nsave::Int

    function SolutionPDAE(d, m, ntime, nsave)
        @assert T <: Real
        @assert d > 0
        @assert m > 0
        @assert nsave > 0
        @assert ntime ≥ nsave
        @assert mod(ntime, nsave) == 0

        n = div(ntime, nsave)
        t = zeros(T, n+1)
        q = zeros(T, d, n+1)
        p = zeros(T, d, n+1)
        λ = zeros(T, m, n+1)
        new(d, m, n, t, q, p, λ, ntime, nsave)
    end
end

function SolutionPDAE(equation::PDAE, ntime::Int, nsave::Int=1)
    T1 = eltype(equation.q₀)
    T2 = eltype(equation.p₀)
    T3 = eltype(equation.λ₀)
    @assert T1 == T2 == T3
    SolutionPDAE{T1}(equation.m, equation.n, ntime, nsave)
end

function Base.getindex(s::SolutionPDAE, i::Int)
    # TODO
end

function reset(s::SolutionPDAE)
    # TODO
end
