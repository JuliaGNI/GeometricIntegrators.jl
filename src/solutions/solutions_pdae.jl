
abstract type SolutionPDAE{dType, tType, N} <: Solution{dType, tType, N} end


"Serial Solution of a partitioned differential algebraic equation."
mutable struct SSolutionPDAE{dType, tType, N} <: SolutionPDAE{dType, tType, N}
    nd::Int
    nm::Int
    nt::Int
    ni::Int
    t::TimeSeries{tType}
    q::SDataSeries{dType,N}
    p::SDataSeries{dType,N}
    λ::SDataSeries{dType,N}
    ntime::Int
    nsave::Int
    counter::Int
end

function SSolutionPDAE(equation::Union{IODE{DT,TT},PDAE{DT,TT},IDAE{DT,TT}}, Δt::TT, ntime::Int, nsave::Int=1) where {DT,TT}
    N  = equation.n > 1 ? 3 : 2
    nd = equation.d
    nm = equation.m
    ni = equation.n
    nt = div(ntime, nsave)

    @assert DT <: Number
    @assert TT <: Real
    @assert nd > 0
    @assert nm > 0
    @assert ni > 0
    @assert nsave > 0
    @assert ntime ≥ nsave
    @assert mod(ntime, nsave) == 0

    t = TimeSeries{TT}(nt, Δt, nsave)
    q = SDataSeries(DT, nd, nt, ni)
    p = SDataSeries(DT, nd, nt, ni)
    λ = SDataSeries(DT, nm, nt, ni)
    s = SSolutionPDAE{DT,TT,N}(nd, nm, nt, ni, t, q, p, λ, ntime, nsave, 0)
    set_initial_conditions!(s, equation)
    return s
end


"Parallel Solution of a partitioned differential algebraic equation."
mutable struct PSolutionPDAE{dType, tType, N} <: SolutionPDAE{dType, tType, N}
    nd::Int
    nm::Int
    nt::Int
    ni::Int
    t::TimeSeries{tType}
    q::PDataSeries{dType,N}
    p::PDataSeries{dType,N}
    λ::PDataSeries{dType,N}
    ntime::Int
    nsave::Int
    counter::Int
end

function PSolutionPDAE(equation::Union{IODE{DT,TT},PDAE{DT,TT},IDAE{DT,TT}}, Δt::TT, ntime::Int, nsave::Int=1) where {DT,TT}
    N  = equation.n > 1 ? 3 : 2
    nd = equation.d
    nm = equation.m
    ni = equation.n
    nt = div(ntime, nsave)

    @assert DT <: Number
    @assert TT <: Real
    @assert nd > 0
    @assert nm > 0
    @assert ni > 0
    @assert nsave > 0
    @assert ntime ≥ nsave
    @assert mod(ntime, nsave) == 0

    t = TimeSeries{TT}(nt, Δt, nsave)
    q = PDataSeries(DT, nd, nt, ni)
    p = PDataSeries(DT, nd, nt, ni)
    λ = PDataSeries(DT, nm, nt, ni)
    s = PSolutionPDAE{DT,TT,N}(nd, nm, nt, ni, t, q, p, λ, ntime, nsave, 0)
    set_initial_conditions!(s, equation)
    return s
end


function set_initial_conditions!(sol::SolutionPDAE{DT,TT}, equ::Union{PDAE{DT,TT},IDAE{DT,TT}}) where {DT,TT}
    set_initial_conditions!(sol, equ.t₀, equ.q₀, equ.p₀, equ.λ₀)
end

function set_initial_conditions!(sol::SolutionPDAE{DT,TT}, equ::IODE{DT,TT}) where {DT,TT}
    set_initial_conditions!(sol, equ.t₀, equ.q₀, equ.p₀, zeros(equ.q₀))
end

function set_initial_conditions!(sol::SolutionPDAE{DT,TT}, t₀::TT, q₀::Union{Array{DT}, Array{Double{DT}}}, p₀::Union{Array{DT}, Array{Double{DT}}}, λ₀::Union{Array{DT}, Array{Double{DT}}}) where {DT,TT}
    set_data!(sol.q, q₀, 0)
    set_data!(sol.p, p₀, 0)
    set_data!(sol.λ, λ₀, 0)
    compute_timeseries!(sol.t, t₀)
end

function get_initial_conditions!(sol::SolutionPDAE{DT,TT}, q::Union{Vector{DT}, Vector{Double{DT}}}, p::Union{Vector{DT}, Vector{Double{DT}}}, λ::Union{Vector{DT}, Vector{Double{DT}}}, k) where {DT,TT}
    get_data!(sol.q, q, 0, k)
    get_data!(sol.p, p, 0, k)
    get_data!(sol.λ, λ, 0, k)
end

function get_initial_conditions!(sol::SolutionPDAE{DT,TT}, q::Union{Vector{DT}, Vector{Double{DT}}}, p::Union{Vector{DT}, Vector{Double{DT}}}, k) where {DT,TT}
    get_data!(sol.q, q, 0, k)
    get_data!(sol.p, p, 0, k)
end

function copy_solution!(sol::SolutionPDAE{DT,TT}, q::Union{Vector{DT}, Vector{Double{DT}}}, p::Union{Vector{DT}, Vector{Double{DT}}}, λ::Union{Vector{DT}, Vector{Double{DT}}}, n, k) where {DT,TT}
    if mod(n, sol.nsave) == 0
        j = div(n, sol.nsave)
        set_data!(sol.q, q, j, k)
        set_data!(sol.p, p, j, k)
        set_data!(sol.λ, λ, j, k)
        sol.counter += 1
    end
end

function copy_solution!(sol::SolutionPDAE{DT,TT}, q::Union{Vector{DT}, Vector{Double{DT}}}, p::Union{Vector{DT}, Vector{Double{DT}}}, n, k) where {DT,TT}
    if mod(n, sol.nsave) == 0
        j = div(n, sol.nsave)
        set_data!(sol.q, q, j, k)
        set_data!(sol.p, p, j, k)
        sol.counter += 1
    end
end

function reset!(s::SolutionPDAE)
    reset!(s.q)
    reset!(s.p)
    reset!(s.λ)
    compute_timeseries!(solution.t, solution.t[end])
    s.counter = 0
end


"Creates HDF5 file and initialises datasets for PDAE solution object."
function create_hdf5(solution::SolutionPDAE{DT,TT}, file::AbstractString, ntime::Int=1) where {DT,TT}
    @assert ntime ≥ 1

    info("Creating HDF5 file ", file)
    isfile(file) ? warn("Overwriting existing HDF5 file.") : nothing
    h5 = h5open(file, "w")

    # save attributes
    attrs(h5)["ntime"] = solution.ntime
    attrs(h5)["nsave"] = solution.nsave

    # copy time
    write(h5, "t", solution.t.t)

    # create datasets
    q = d_create(h5, "q", datatype(DT), dataspace(solution.nd, solution.nt+1), "chunk", (solution.nd,1))
    p = d_create(h5, "p", datatype(DT), dataspace(solution.nd, solution.nt+1), "chunk", (solution.nd,1))
    λ = d_create(h5, "λ", datatype(DT), dataspace(solution.nm, solution.nt+1), "chunk", (solution.nm,1))

    # copy initial conditions
    q[1:solution.nd, 1] = solution.q.d[1:solution.nd, 1]
    p[1:solution.nd, 1] = solution.p.d[1:solution.nd, 1]
    λ[1:solution.nm, 1] = solution.λ.d[1:solution.nm, 1]

    return h5
end

"Append solution to HDF5 file."
function CommonFunctions.write_to_hdf5(solution::SolutionPDAE, h5::HDF5.HDF5File, offset=0)
    # aquire dataset from HDF5 file
    q = h5["q"]
    p = h5["p"]
    λ = h5["λ"]

    # set convenience variables and compute ranges
    d  = solution.nd
    m  = solution.nm
    n  = solution.nt
    j1 = offset+2
    j2 = offset+1+n

    # copy data from solution to HDF5 dataset
    q[1:d, j1:j2] = solution.q.d[1:d, 2:n+1]
    p[1:d, j1:j2] = solution.p.d[1:d, 2:n+1]
    λ[1:m, j1:j2] = solution.λ.d[1:m, 2:n+1]

    return nothing
end
