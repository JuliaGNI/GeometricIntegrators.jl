
"Solution of a partitioned ordinary differential equation."
struct SolutionPODE{dType, tType, N} <: Solution{dType, tType, N}
    nd::Int
    nt::Int
    ni::Int
    t::TimeSeries{tType}
    q::SDataSeries{dType,N}
    p::SDataSeries{dType,N}
    ntime::Int
    nsave::Int
end

function SolutionPODE(equation::Union{PODE{DT,TT}, IODE{DT,TT}}, Δt::TT, ntime::Int, nsave::Int=1) where {DT,TT}
    N  = equation.n > 1 ? 3 : 2
    nd = equation.d
    ni = equation.n
    nt = div(ntime, nsave)

    @assert DT <: Number
    @assert TT <: Real
    @assert nd > 0
    @assert ni > 0
    @assert nsave > 0
    @assert ntime ≥ nsave
    @assert mod(ntime, nsave) == 0

    t = TimeSeries{TT}(nt, Δt, nsave)
    q = SDataSeries(DT, nd, nt, ni)
    p = SDataSeries(DT, nd, nt, ni)
    s = SolutionPODE{DT,TT,N}(nd, nt, ni, t, q, p, ntime, nsave)
    set_initial_conditions!(s, equation)
    return s
end

function set_initial_conditions!(sol::SolutionPODE{DT,TT}, equ::Union{PODE{DT,TT},IODE{DT,TT}}) where {DT,TT}
    set_initial_conditions!(sol, equ.t₀, equ.q₀, equ.p₀)
end

function set_initial_conditions!(sol::SolutionPODE{DT,TT}, t₀::TT, q₀::Array{DT}, p₀::Array{DT}) where {DT,TT}
    set_data!(sol.q, q₀, 0)
    set_data!(sol.p, p₀, 0)
    compute_timeseries!(sol.t, t₀)
end

function get_initial_conditions!(sol::SolutionPODE{DT,TT}, q::Vector{DT}, p::Vector{DT}, k) where {DT,TT}
    get_data!(sol.q, q, 0, k)
    get_data!(sol.p, p, 0, k)
end

function copy_solution!(sol::SolutionPODE{DT,TT}, q::Vector{DT}, p::Vector{DT}, n, k) where {DT,TT}
    if mod(n, sol.nsave) == 0
        j = div(n, sol.nsave)
        set_data!(sol.q, q, j, k)
        set_data!(sol.p, p, j, k)
    end
end

function reset!(s::SolutionPODE)
    reset!(s.q)
    reset!(s.p)
    compute_timeseries!(solution.t, solution.t[end])
end
