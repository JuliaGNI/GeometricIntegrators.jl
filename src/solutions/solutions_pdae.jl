
"Solution of a partitioned differential algebraic equation."
immutable SolutionPDAE{dType, tType, N} <: Solution{dType, tType, N}
    nd::Int
    nm::Int
    nt::Int
    n0::Int
    t::TimeSeries{tType}
    x::Array{dType, N}
    q::AbstractArray{dType} # TODO Provide actual type.
    p::AbstractArray{dType} # TODO Provide actual type.
    λ::AbstractArray{dType} # TODO Provide actual type.
    ntime::Int
    nsave::Int

    function SolutionPDAE(nd, nm, n0, ntime, nsave, Δt)
        @assert dType <: Number
        @assert tType <: Real
        @assert nd > 0
        @assert nm > 0
        @assert n0 > 0
        @assert nsave > 0
        @assert ntime ≥ nsave
        @assert mod(ntime, nsave) == 0

        nt = div(ntime, nsave)
        t = TimeSeries{tType}(nt, Δt, nsave)

        @assert N ∈ (3,4)

        if N == 3
            x = zeros(dType, 3, nd, nt+1)
            q = view(x, 1, :, :)
            p = view(x, 2, :, :)
            λ = view(x, 3, :, :)
        elseif N == 4
            x = zeros(dType, 3, nd, nt+1, n0)
            q = view(x, 1, :, :, :)
            p = view(x, 2, :, :, :)
            λ = view(x, 3, :, :, :)
        end

        new(nd, nm, nt, n0, t, x, q, p, λ, ntime, nsave)
    end
end

function SolutionPDAE{DT,TT}(equation::PDAE{DT,TT}, Δt::TT, ntime::Int, nsave::Int=1)
    N = equation.n > 1 ? 4 : 3
    s = SolutionPDAE{DT,TT,N}(equation.d, equation.m, equation.n, ntime, nsave, Δt)
    set_initial_conditions!(s, equation)
    return s
end

function set_initial_conditions!{DT,TT}(solution::SolutionPDAE{DT,TT}, equation::PDAE{DT,TT})
    set_initial_conditions!(solution, equation.t₀, equation.q₀, equation.p₀, equation.λ₀)
end

function set_initial_conditions!{DT,TT}(solution::SolutionPDAE{DT,TT,3}, t₀::TT, q₀::Array{DT,1}, p₀::Array{DT,1}, λ₀::Array{DT,1})
    @assert size(solution, 2) == size(q₀, 1) == size(p₀, 1) ≥ size(λ₀, 1)
    @inbounds for i in 1:size(solution, 2)
        solution[1, i, 0] = q₀[i]
        solution[2, i, 0] = p₀[i]
    end
    for i in 1:size(λ₀,1)
        solution[3, i, 0] = λ₀[i]
    end
    compute_timeseries!(solution.t, t₀)
end

function set_initial_conditions!{DT,TT}(solution::SolutionPDAE{DT,TT,4}, t₀::TT, q₀::Array{DT,2}, p₀::Array{DT,2}, λ₀::Array{DT,2})
    @assert size(solution, 2) == size(q₀, 1) == size(p₀, 1) ≥ size(λ₀, 1)
    @assert size(solution, 4) == size(q₀, 2) == size(p₀, 2) ≥ size(λ₀, 2)
    @inbounds for k in 1:size(solution,4)
        for i in 1:size(solution,2)
            solution[1, i, 0, k] = q₀[i,k]
            solution[2, i, 0, k] = p₀[i,k]
        end
        for i in 1:size(λ₀,1)
            solution[3, i, 0, k] = λ₀[i,k]
        end
    end
    compute_timeseries!(solution.t, t₀)
end


function copy_solution!{DT,TT}(q::Vector{DT}, p::Vector{DT}, λ::Vector{DT}, sol::SolutionPDAE{DT,TT,3}, n, k)
    if mod(n, sol.nsave) == 0
        j = div(n, sol.nsave)+1
        @assert size(sol.x, 2) == length(q) == length(p) ≥ length(λ)
        @assert j ≤ size(sol.x, 3)
        @assert k == 1
        @inbounds for i in 1:size(sol.x, 2)
            sol.x[1, i, j] = q[i]
            sol.x[2, i, j] = p[i]
            sol.x[3, i, j] = λ[i]
        end
    end
end

function copy_solution!{DT,TT}(q::Vector{DT}, p::Vector{DT}, λ::Vector{DT}, sol::SolutionPDAE{DT,TT,4}, n, k)
    if mod(n, sol.nsave) == 0
        j = div(n, sol.nsave)+1
        @assert size(sol.x, 2) == length(q) == length(p) ≥ length(λ)
        @assert j ≤ size(sol.x, 3)
        @assert k ≤ size(sol.x, 4)
        @inbounds for i in 1:size(sol.x, 2)
            sol.x[1, i, j, k] = q[i]
            sol.x[2, i, j, k] = p[i]
            sol.x[3, i, j, k] = λ[i]
        end
    end
end

function reset!{DT,TT}(s::SolutionPDAE{DT,TT,3})
    for j in 1:size(solution,2)
        for i in 1:size(solution,1)
            solution[i, j, 0] = solution[i, j, end]
        end
    end
end

function reset!{DT,TT}(s::SolutionPDAE{DT,TT,4})
    for k in 1:size(solution,4)
        for j in 1:size(solution,2)
            for i in 1:size(solution,1)
                solution[i, j, 0, k] = solution[i, j, end, k]
            end
        end
    end
end

Base.indices{DT,TT}(s::SolutionPDAE{DT,TT,3}) = (1:3, 1:s.nd, 0:s.nt)
Base.indices{DT,TT}(s::SolutionPDAE{DT,TT,4}) = (1:3, 1:s.nd, 0:s.nt, 1:s.n0)
Base.strides(s::SolutionPDAE) = strides(s.x)

@inline function Base.getindex{DT,TT}(s::SolutionPDAE{DT,TT,3}, j::Int, k::Int, m::Int)
    @boundscheck checkbounds(s.x, j, k, m+1)
    @inbounds r = getindex(s.x, j, k, m+1)
    return r
end

@inline function Base.getindex{DT,TT}(s::SolutionPDAE{DT,TT,3}, k::Int, m::Int)
    @boundscheck checkbounds(s.x, :, k, m+1)
    @inbounds r = getindex(s.x, :, k, m+1)
    return r
end

@inline function Base.getindex{DT,TT}(s::SolutionPDAE{DT,TT,3}, m::Int)
    @boundscheck checkbounds(s.x, :, :, m)
    @inbounds r = getindex(s.x, :, :, m)
    return r
end

@inline function Base.getindex{DT,TT}(s::SolutionPDAE{DT,TT,4}, i::Int, j::Int, k::Int, m::Int)
    @boundscheck checkbounds(s.x, i, j, k+1, m)
    @inbounds r = getindex(s.x, i, j, k+1, m)
    return r
end

@inline function Base.getindex{DT,TT}(s::SolutionPDAE{DT,TT,4}, j::Int, k::Int, m::Int)
    @boundscheck checkbounds(s.x, :, j, k+1, m)
    @inbounds r = getindex(s.x, :, j, k+1, m)
    return r
end

@inline function Base.getindex{DT,TT}(s::SolutionPDAE{DT,TT,4}, k::Int, m::Int)
    @boundscheck checkbounds(s.x, :, :, k+1, m)
    @inbounds r = getindex(s.x, :, :, k+1, m)
    return r
end

@inline function Base.getindex{DT,TT}(s::SolutionPDAE{DT,TT,4}, m::Int)
    @boundscheck checkbounds(s.x, :, :, :, m)
    @inbounds r = getindex(s.x, :, :, :, m)
    return r
end

@inline function Base.setindex!{DT,TT}(s::SolutionPDAE{DT,TT,3}, x, j::Int, k::Int, m::Int)
    @assert length(x) == 1
    @boundscheck checkbounds(s.x, j, k, m+1)
    @inbounds setindex!(s.x, x, j, k, m+1)
end

@inline function Base.setindex!{DT,TT}(s::SolutionPDAE{DT,TT,3}, x, k::Int, m::Int)
    @assert ndims(x) == 1
    @assert length(x) == 3
    @boundscheck checkbounds(s.x, :, k, m+1)
    @inbounds setindex!(s.x, x, :, k, m+1)
end

@inline function Base.setindex!{DT,TT}(s::SolutionPDAE{DT,TT,3}, x, m::Int)
    @assert ndims(x) == 2
    @assert size(x, 1) == size(s.x, 1)
    @assert size(x, 2) == size(s.x, 2)
    @boundscheck checkbounds(s.x, :, :, m+1)
    @inbounds setindex!(s.x, x, :, :, m+1)
end
@inline function Base.setindex!{DT,TT}(s::SolutionPDAE{DT,TT,4}, x, i::Int, j::Int, k::Int, m::Int)
    @assert length(x) == 1
    @boundscheck checkbounds(s.x, i, j, k+1, m)
    @inbounds setindex!(s.x, x, i, j, k+1, m)
end

@inline function Base.setindex!{DT,TT}(s::SolutionPDAE{DT,TT,4}, x, j::Int, k::Int, m::Int)
    @assert ndims(x) == 1
    @assert length(x) == 3
    @boundscheck checkbounds(s.x, :, j, k+1, m)
    @inbounds setindex!(s.x, x, :, j, k+1, m)
end

@inline function Base.setindex!{DT,TT}(s::SolutionPDAE{DT,TT,4}, x, k::Int, m::Int)
    @assert ndims(x) == 2
    @assert size(x, 1) == size(s.x, 1)
    @assert size(x, 2) == size(s.x, 2)
    @boundscheck checkbounds(s.x, :, :, k+1, m)
    @inbounds setindex!(s.x, x, :, :, k+1, m)
end

@inline function Base.setindex!{DT,TT}(s::SolutionPDAE{DT,TT,4}, x, m::Int)
    @assert ndims(x) == 3
    @assert size(x, 1) == size(s.x, 1)
    @assert size(x, 2) == size(s.x, 2)
    @assert size(x, 3) == size(s.x, 3)
    @boundscheck checkbounds(s.x, :, :, :, m)
    @inbounds setindex!(s.x, x, :, :, :, m)
end
