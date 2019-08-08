
"Solution of a partitioned ordinary differential equation."
mutable struct SolutionPODE{dType, tType, N} <: DeterministicSolution{dType, tType, N}
    nd::Int
    nt::Int
    ni::Int
    t::TimeSeries{tType}
    q::SDataSeries{dType,N}
    p::SDataSeries{dType,N}
    ntime::Int
    nsave::Int
    counter::Int
end

function SolutionPODE(equation::Union{PODE{DT,TT}, IODE{DT,TT}, VODE{DT,TT}}, Δt::TT, ntime::Int, nsave::Int=1) where {DT,TT}
    N  = equation.n > 1 ? 3 : 2
    nd = equation.d
    ni = equation.n
    nt = div(ntime, nsave)

    @assert DT <: Number
    @assert TT <: Real
    @assert nd > 0
    @assert ni > 0
    @assert nsave > 0
    @assert ntime == 0 || ntime ≥ nsave
    @assert mod(ntime, nsave) == 0

    t = TimeSeries{TT}(nt, Δt, nsave)
    q = SDataSeries(DT, nd, nt, ni)
    p = SDataSeries(DT, nd, nt, ni)
    s = SolutionPODE{DT,TT,N}(nd, nt, ni, t, q, p, ntime, nsave, 0)
    set_initial_conditions!(s, equation)
    return s
end

function SolutionPODE(t::TimeSeries{TT}, q::SDataSeries{DT,N}, p::SDataSeries{DT,N}, ntime::Int, nsave::Int) where {DT,TT,N}
    @assert q.nd == p.nd
    @assert q.nt == p.nt
    @assert q.ni == p.ni

    # extract parameters
    nd = q.nd
    ni = q.ni
    nt = t.n

    # create solution
    SolutionPODE{DT,TT,N}(nd, nt, ni, t, q, p, ntime, nsave, 0)
end

function SolutionPODE(file::String)
    # open HDF5 file
    info("Reading HDF5 file ", file)
    h5 = h5open(file, "r")

    # read attributes
    ntime = read(attrs(h5)["ntime"])
    nsave = read(attrs(h5)["nsave"])

    # reading data arrays
    t = TimeSeries(read(h5["t"]), nsave)
    q = SDataSeries(read(h5["q"]))
    p = SDataSeries(read(h5["p"]))

    # create solution
    SolutionPODE(t, q, p, ntime, nsave)
end

Base.:(==)(sol1::SolutionPODE, sol2::SolutionPODE) = (
                                sol1.nd == sol2.nd
                             && sol1.nt == sol2.nt
                             && sol1.ni == sol2.ni
                             && sol1.t  == sol2.t
                             && sol1.q  == sol2.q
                             && sol1.p  == sol2.p
                             && sol1.ntime == sol2.ntime
                             && sol1.nsave == sol2.nsave
                             && sol1.counter == sol2.counter)

time(sol::SolutionPODE)  = sol.t.t
ntime(sol::SolutionPODE) = sol.ntime
nsave(sol::SolutionPODE) = sol.nsave


function set_initial_conditions!(sol::SolutionPODE{DT,TT}, equ::Union{PODE{DT,TT},IODE{DT,TT},VODE{DT,TT}}) where {DT,TT}
    set_initial_conditions!(sol, equ.t₀, equ.q₀, equ.p₀)
end

function set_initial_conditions!(sol::SolutionPODE{DT,TT}, t₀::TT, q₀::Union{Array{DT}, Array{TwicePrecision{DT}}}, p₀::Union{Array{DT}, Array{TwicePrecision{DT}}}) where {DT,TT}
    set_data!(sol.q, q₀, 0)
    set_data!(sol.p, p₀, 0)
    compute_timeseries!(sol.t, t₀)
end

function get_initial_conditions!(sol::SolutionPODE{DT,TT}, q::Union{Vector{DT}, Vector{TwicePrecision{DT}}}, p::Union{Vector{DT}, Vector{TwicePrecision{DT}}}, k) where {DT,TT}
    get_data!(sol.q, q, 0, k)
    get_data!(sol.p, p, 0, k)
end

function copy_solution!(sol::SolutionPODE{DT,TT}, q::Union{Vector{DT}, Vector{TwicePrecision{DT}}}, p::Union{Vector{DT}, Vector{TwicePrecision{DT}}}, n, k) where {DT,TT}
    if mod(n, sol.nsave) == 0
        j = div(n, sol.nsave)
        set_data!(sol.q, q, j, k)
        set_data!(sol.p, p, j, k)
        sol.counter += 1
    end
end

function reset!(sol::SolutionPODE)
    reset!(sol.q)
    reset!(sol.p)
    compute_timeseries!(sol.t, sol.t[end])
    sol.counter = 0
end


"Creates HDF5 file and initialises datasets for PODE solution object."
function create_hdf5(solution::SolutionPODE{DT,TT,2}, file::AbstractString, ntime::Int=1) where {DT,TT}
    @assert ntime ≥ 1

    # create HDF5 file and save attributes and common parameters
    h5 = createHDF5(solution, file)

    # create datasets
    q = d_create(h5, "q", datatype(DT), dataspace(solution.nd, solution.nt+1), "chunk", (solution.nd,1))
    p = d_create(h5, "p", datatype(DT), dataspace(solution.nd, solution.nt+1), "chunk", (solution.nd,1))

    # copy initial conditions
    q[1:solution.nd, 1] = solution.q.d[1:solution.nd, 1]
    p[1:solution.nd, 1] = solution.p.d[1:solution.nd, 1]

    return h5
end

"Creates HDF5 file and initialises datasets for PODE solution object."
function create_hdf5(solution::SolutionPODE{DT,TT,3}, file::AbstractString, ntime::Int=1) where {DT,TT}
    @assert ntime ≥ 1

    # create HDF5 file and save attributes and common parameters
    h5 = createHDF5(solution, file)

    # create datasets
    q = d_create(h5, "q", datatype(DT), dataspace(solution.nd, solution.nt+1, solution.ni), "chunk", (solution.nd,1,1))
    p = d_create(h5, "p", datatype(DT), dataspace(solution.nd, solution.nt+1, solution.ni), "chunk", (solution.nd,1,1))

    # copy initial conditions
    q[1:solution.nd, 1, 1:solution.ni] = solution.q.d[1:solution.nd, 1, 1:solution.ni]
    p[1:solution.nd, 1, 1:solution.ni] = solution.p.d[1:solution.nd, 1, 1:solution.ni]

    return h5
end

"Append solution to HDF5 file."
function CommonFunctions.write_to_hdf5(solution::SolutionPODE{DT,TT,2}, h5::HDF5.HDF5File, offset=0) where {DT,TT}
    # set convenience variables and compute ranges
    d  = solution.nd
    n  = solution.nt
    j1 = offset+2
    j2 = offset+1+n

    # copy data from solution to HDF5 dataset
    h5["q"][1:d, j1:j2] = solution.q.d[1:d, 2:n+1]
    h5["p"][1:d, j1:j2] = solution.p.d[1:d, 2:n+1]

    return nothing
end

"Append solution to HDF5 file."
function CommonFunctions.write_to_hdf5(solution::SolutionPODE{DT,TT,3}, h5::HDF5.HDF5File, offset=0) where {DT,TT}
    # set convenience variables and compute ranges
    d  = solution.nd
    n  = solution.nt
    i  = solution.ni
    j1 = offset+2
    j2 = offset+1+n

    # copy data from solution to HDF5 dataset
    h5["q"][1:d, j1:j2, 1:i] = solution.q.d[1:d, 2:n+1, 1:i]
    h5["p"][1:d, j1:j2, 1:i] = solution.p.d[1:d, 2:n+1, 1:i]

    return nothing
end
