
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

function set_initial_conditions!(sol::SolutionPODE{DT,TT}, t₀::TT, q₀::Union{Array{DT}, Array{Double{DT}}}, p₀::Union{Array{DT}, Array{Double{DT}}}) where {DT,TT}
    set_data!(sol.q, q₀, 0)
    set_data!(sol.p, p₀, 0)
    compute_timeseries!(sol.t, t₀)
end

function get_initial_conditions!(sol::SolutionPODE{DT,TT}, q::Union{Vector{DT}, Vector{Double{DT}}}, p::Union{Vector{DT}, Vector{Double{DT}}}, k) where {DT,TT}
    get_data!(sol.q, q, 0, k)
    get_data!(sol.p, p, 0, k)
end

function copy_solution!(sol::SolutionPODE{DT,TT}, q::Union{Vector{DT}, Vector{Double{DT}}}, p::Union{Vector{DT}, Vector{Double{DT}}}, n, k) where {DT,TT}
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


"Creates HDF5 file and initialises datasets for PODE solution object."
function create_hdf5(solution::SolutionPODE{DT,TT}, file::AbstractString, ntime::Int=1) where {DT,TT}
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

    # copy initial conditions
    q[1:solution.nd, 1] = solution.q.d[1:solution.nd, 1]
    p[1:solution.nd, 1] = solution.p.d[1:solution.nd, 1]

    return h5
end

"Append solution to HDF5 file."
function CommonFunctions.write_to_hdf5(solution::SolutionPODE, h5::HDF5.HDF5File, offset=0)
    # aquire dataset from HDF5 file
    t = h5["t"]
    q = h5["q"]
    p = h5["p"]

    # set convenience variables and compute ranges
    d  = solution.nd
    n  = solution.nt
    j1 = offset+2
    j2 = offset+1+n

    # copy data from solution to HDF5 dataset
    q[1:d, j1:j2] = solution.q.d[1:d, 2:n+1]
    p[1:d, j1:j2] = solution.p.d[1:d, 2:n+1]

    return nothing
end
