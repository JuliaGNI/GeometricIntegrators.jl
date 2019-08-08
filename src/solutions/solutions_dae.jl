
"Solution of a differential algebraic equation."
struct SolutionDAE{dType, tType, N} <: DeterministicSolution{dType, tType, N}
    nd::Int
    nm::Int
    nt::Int
    ni::Int
    t::TimeSeries{tType}
    q::SDataSeries{dType,N}
    λ::SDataSeries{dType,N}
    ntime::Int
    nsave::Int
end

function SolutionDAE(equation::DAE{DT,TT}, Δt::TT, ntime::Int, nsave::Int=1) where {DT,TT}
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
    @assert ntime == 0 || ntime ≥ nsave
    @assert mod(ntime, nsave) == 0

    t = TimeSeries{TT}(nt, Δt, nsave)
    q = SDataSeries(DT, nd, nt, ni)
    λ = SDataSeries(DT, nm, nt, ni)
    s = SolutionDAE{DT,TT,N}(nd, nm, nt, ni, t, q, λ, ntime, nsave)
    set_initial_conditions!(s, equation)
    return s
end

function SolutionDAE(t::TimeSeries{TT}, q::SDataSeries{DT,N}, λ::SDataSeries{DT,N}, ntime::Int, nsave::Int) where {DT,TT,N}
    @assert q.nd
    @assert q.nt == λ.nt
    @assert q.ni == λ.ni

    # extract parameters
    nd = q.nd
    nm = λ.nd
    ni = q.ni
    nt = t.n

    # create solution
    SolutionDAE{DT,TT,N}(nd, nm, nt, ni, t, q, λ, ntime, nsave, 0)
end

function SolutionDAE(file::String)
    # open HDF5 file
    info("Reading HDF5 file ", file)
    h5 = h5open(file, "r")

    # read attributes
    ntime = read(attrs(h5)["ntime"])
    nsave = read(attrs(h5)["nsave"])

    # reading data arrays
    t = TimeSeries(read(h5["t"]), nsave)
    q = SDataSeries(read(h5["q"]))
    λ = SDataSeries(read(h5["λ"]))

    # create solution
    SolutionDAE(t, q, λ, ntime, nsave)
end


time(sol::SolutionDAE)  = sol.t.t
ntime(sol::SolutionDAE) = sol.ntime
nsave(sol::SolutionDAE) = sol.nsave


function set_initial_conditions!(sol::SolutionDAE{DT,TT}, equ::DAE{DT,TT}) where {DT,TT}
    set_initial_conditions!(sol, equ.t₀, equ.q₀, equ.λ₀)
end

function set_initial_conditions!(sol::SolutionDAE{DT,TT}, t₀::TT, q₀::Array{DT}, λ₀::Array{DT}) where {DT,TT}
    set_data!(sol.q, q₀, 0)
    set_data!(sol.λ, λ₀, 0)
    compute_timeseries!(sol.t, t₀)
end

function get_initial_conditions!(sol::SolutionDAE{DT,TT}, q::Vector{DT}, λ::Vector{DT}, k) where {DT,TT}
    get_data!(sol.q, q, 0, k)
    get_data!(sol.λ, λ, 0, k)
end

function copy_solution!(sol::SolutionDAE{DT,TT}, q::Vector{DT}, λ::Vector{DT}, n, k) where {DT,TT}
    if mod(n, sol.nsave) == 0
        j = div(n, sol.nsave)
        set_data!(sol.q, q, j, k)
        set_data!(sol.λ, λ, j, k)
    end
end

function reset!(sol::SolutionDAE)
    reset!(sol.q)
    reset!(sol.λ)
    compute_timeseries!(sol.t, sol.t[end])
end


"Creates HDF5 file and initialises datasets for DAE solution object."
function create_hdf5(solution::SolutionDAE{DT,TT,2}, file::AbstractString, ntime::Int=1) where {DT,TT}
    @assert ntime ≥ 1

    # create HDF5 file and save attributes and common parameters
    h5 = createHDF5(solution, file)

    # save attributes
    save_attributes(solution, h5)

    # create datasets
    q = d_create(h5, "q", datatype(DT), dataspace(solution.nd, solution.nt+1), "chunk", (solution.nd,1))
    λ = d_create(h5, "λ", datatype(DT), dataspace(solution.nm, solution.nt+1), "chunk", (solution.nm,1))

    # copy initial conditions
    q[1:solution.nd, 1] = solution.q.d[1:solution.nd, 1]
    λ[1:solution.nm, 1] = solution.λ.d[1:solution.nm, 1]

    return h5
end

"Creates HDF5 file and initialises datasets for DAE solution object."
function create_hdf5(solution::SolutionDAE{DT,TT,3}, file::AbstractString, ntime::Int=1) where {DT,TT}
    @assert ntime ≥ 1

    # create HDF5 file and save attributes and common parameters
    h5 = createHDF5(solution, file)

    # save attributes
    save_attributes(solution, h5)

    # create datasets
    q = d_create(h5, "q", datatype(DT), dataspace(solution.nd, solution.nt+1, solution.ni), "chunk", (solution.nd,1,1))
    λ = d_create(h5, "λ", datatype(DT), dataspace(solution.nm, solution.nt+1, solution.ni), "chunk", (solution.nm,1,1))

    # copy initial conditions
    q[1:solution.nd, 1, 1:solution.ni] = solution.q.d[1:solution.nd, 1, 1:solution.ni]
    λ[1:solution.nm, 1, 1:solution.ni] = solution.λ.d[1:solution.nm, 1, 1:solution.ni]

    return h5
end

"Append solution to HDF5 file."
function CommonFunctions.write_to_hdf5(solution::SolutionDAE{DT,TT,2}, h5::HDF5.HDF5File, offset=0) where {DT,TT}
    # set convenience variables and compute ranges
    d  = solution.nd
    m  = solution.nm
    n  = solution.nt
    j1 = offset+2
    j2 = offset+1+n

    # copy data from solution to HDF5 dataset
    h5["q"][1:d, j1:j2] = solution.q.d[1:d, 2:n+1]
    h5["λ"][1:m, j1:j2] = solution.λ.d[1:m, 2:n+1]

    return nothing
end

"Append solution to HDF5 file."
function CommonFunctions.write_to_hdf5(solution::SolutionDAE{DT,TT,3}, h5::HDF5.HDF5File, offset=0) where {DT,TT}
    # set convenience variables and compute ranges
    d  = solution.nd
    m  = solution.nm
    n  = solution.nt
    i  = solution.ni
    j1 = offset+2
    j2 = offset+1+n

    # copy data from solution to HDF5 dataset
    h5["q"][1:d, j1:j2, 1:i] = solution.q.d[1:d, 2:n+1, 1:i]
    h5["λ"][1:m, j1:j2, 1:i] = solution.λ.d[1:m, 2:n+1, 1:i]

    return nothing
end
