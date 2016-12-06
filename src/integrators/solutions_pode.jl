
"Solution of a partitioned ordinary differential equation."
immutable SolutionPODE{dType, tType} <: Solution{dType, tType, 4}
    nd::Int
    nt::Int
    n0::Int
    t::Timeseries{tType}
    x::Array{dType, 4}
    q::AbstractArray{dType,3}
    p::AbstractArray{dType,3}
    ntime::Int
    nsave::Int

    function SolutionPODE(nd, n0, ntime, nsave, Δt)
        @assert dType <: Number
        @assert tType <: Real
        @assert nd > 0
        @assert n0 > 0
        @assert nsave > 0
        @assert ntime ≥ nsave
        @assert mod(ntime, nsave) == 0

        nt = div(ntime, nsave)
        t = Timeseries{tType}(nt, Δt, nsave)
        x = zeros(dType, 2, nd, nt+1, n0)
        q = view(x, 1, :, :, :)
        p = view(x, 2, :, :, :)
        new(nd, nt, n0, t, x, q, p, ntime, nsave)
    end
end

function SolutionPODE{DT,TT}(equation::Union{PODE{DT,TT}, IODE{DT}}, Δt::TT, ntime::Int, nsave::Int=1)
    s = SolutionPODE{DT,TT}(equation.d, equation.n, ntime, nsave, Δt)
    set_initial_conditions!(s, equation)
    return s
end

function set_initial_conditions!(solution::SolutionPODE, equation::Union{PODE,IODE})
    for k in 1:solution.n0
        for i in 1:solution.nd
            solution[1, i, 0, k] = equation.q₀[i,k]
            solution[2, i, 0, k] = equation.p₀[i,k]
        end
    end
    solution.t[0] = equation.t₀
    compute_timeseries!(solution.t)
end

function reset!(s::SolutionPODE)
    for k in 1:solution.n0
        for i in 1:solution.nd
            solution[1, i, 0, k] = solution[1, i, end, k]
            solution[2, i, 0, k] = solution[2, i, end, k]
        end
    end
end

Base.indices(s::SolutionPODE) = (1:2, 1:s.nd, 0:s.nt, 1:s.n0)
Base.strides(s::SolutionPODE) = (1, 2, 2*s.nd, 2*s.nd*s.nt)
# Base.linearindexing{T<:SolutionPODE}(::Type{T}) = LinearFast()

@inline function Base.getindex(s::SolutionPODE, i::Int, j::Int, k::Int, m::Int)
    @boundscheck checkbounds(s.x, i, j, k+1, m)
    @inbounds r = getindex(s.x, i, j, k+1, m)
    return r
end

@inline function Base.getindex(s::SolutionPODE, j::Int, k::Int, m::Int)
    if s.n0 == 1
        @boundscheck checkbounds(s.x, j, k, m+1, 1)
        @inbounds r = getindex(s.x, j, k, m+1, 1)
    else
        @boundscheck checkbounds(s.x, :, j, k+1, m)
        @inbounds r = getindex(s.x, :, j, k+1, m)
    end
    return r
end

@inline function Base.getindex(s::SolutionPODE, k::Int, m::Int)
    if s.n0 == 1
        @boundscheck checkbounds(s.x, :, k, m+1, 1)
        @inbounds r = getindex(s.x, :, k, m+1, 1)
    else
        @boundscheck checkbounds(s.x, :, 1:s.nd, k+1, m)
        @inbounds r = getindex(s.x, :, 1:s.nd, k+1, m)
    end
    return r
end

@inline function Base.getindex(s::SolutionPODE, m::Int)
    if s.n0 == 1
        @boundscheck checkbounds(s.x, :, 1:s.nd, m, 1)
        @inbounds r = getindex(s.x, :, 1:s.nd, m, 1)
    else
        @boundscheck checkbounds(s.x, :, 1:s.nd, :, m)
        @inbounds r = getindex(s.x, :, 1:s.nd, :, m)
    end
    return r
end

@inline function Base.setindex!(s::SolutionPODE, x, i::Int, j::Int, k::Int, m::Int)
    @boundscheck checkbounds(s.x, i, j, k+1, m)
    @inbounds setindex!(s.x, x, i, j, k+1, m)
end

@inline function Base.setindex!(s::SolutionPODE, x, j::Int, k::Int, m::Int)
    if s.n0 == 1
        @assert length(x) == 1
        @boundscheck checkbounds(s.x, j, k, m+1, 1)
        @inbounds setindex!(s.x, x, j, k, m+1, 1)
    else
        @assert length(x) == 2
        @boundscheck checkbounds(s.x, :, j, k+1, m)
        @inbounds setindex!(s.x, x, :, j, k+1, m)
    end
end

@inline function Base.setindex!(s::SolutionPODE, x, k::Int, m::Int)
    if s.n0 == 1
        @assert length(x) == 2
        @boundscheck checkbounds(s.x, :, k, m+1, 1)
        @inbounds setindex!(s.x, x, :, k, m+1, 1)
    else
        @assert size(x, 1) == 2
        @assert size(x, 2) == s.nd
        @boundscheck checkbounds(s.x, :, :, k+1, m)
        @inbounds setindex!(s.x, x, :, :, k+1, m)
    end
end

@inline function Base.setindex!(s::SolutionPODE, x, m::Int)
    if s.n0 == 1
        @assert size(x, 1) == 2
        @assert size(x, 2) == s.nd
        @boundscheck checkbounds(s.x, :, :, m+1, 1)
        @inbounds setindex!(s.x, x, :, :, m+1, 1)
    else
        @assert size(x, 1) == 2
        @assert size(x, 2) == s.nd
        @assert size(x, 3) == s.nt
        @boundscheck checkbounds(s.x, :, :, :, m)
        @inbounds setindex!(s.x, x, :, :, :, m)
    end
end
