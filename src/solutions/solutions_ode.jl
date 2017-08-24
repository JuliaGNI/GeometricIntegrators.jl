"""
`SolutionODE`: Solution of an ordinary differential equation

Contains all fields necessary to store the solution of an ODE.

### Fields

* `nd`: dimension of the dynamical variable ``q``
* `nt`: number of time steps to store
* `ni`: number of initial conditions
* `t`:  time steps
* `q`:  solution `q[nd, nt+1, ni]` with `q[:,0,:]` the initial conditions
* `ntime`: number of time steps to compute
* `nsave`: save every nsave'th time step

"""
mutable struct SolutionODE{dType, tType, N} <: Solution{dType, tType, N}
    nd::Int
    nt::Int
    ni::Int
    t::TimeSeries{tType}
    q::SDataSeries{dType,N}
    ntime::Int
    nsave::Int
    counter::Int
end

function SolutionODE(equation::Union{ODE{DT,TT,FT},SODE{DT,TT,FT}}, Δt::TT, ntime::Int, nsave::Int=1) where {DT,TT,FT}
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
    s = SolutionODE{DT,TT,N}(nd, nt, ni, t, q, ntime, nsave, 0)
    set_initial_conditions!(s, equation)
    return s
end

function SolutionODE(t::TimeSeries{TT}, q::SDataSeries{DT,N}, ntime::Int, nsave::Int) where {DT,TT,N}
    # extract parameters
    nd = q.nd
    ni = q.ni
    nt = t.n

    # create solution
    SolutionODE{DT,TT,N}(nd, nt, ni, t, q, ntime, nsave, 0)
end

function SolutionODE(file::String)
    # open HDF5 file
    info("Reading HDF5 file ", file)
    h5 = h5open(file, "r")

    # read attributes
    ntime = read(attrs(h5)["ntime"])
    nsave = read(attrs(h5)["nsave"])

    # reading data arrays
    t = TimeSeries(read(h5["t"]), nsave)
    q = SDataSeries(read(h5["q"]))

    # create solution
    SolutionODE(t, q, ntime, nsave)
end


time(sol::SolutionODE)  = sol.t.t
ntime(sol::SolutionODE) = sol.ntime
nsave(sol::SolutionODE) = sol.nsave


function set_initial_conditions!(sol::SolutionODE{DT,TT}, equ::Union{ODE{DT,TT},SODE{DT,TT}}) where {DT,TT}
    set_initial_conditions!(sol, equ.t₀, equ.q₀)
end

function set_initial_conditions!(sol::SolutionODE{DT,TT}, t₀::TT, q₀::Union{Array{DT}, Array{Double{DT}}}) where {DT,TT}
    set_data!(sol.q, q₀, 0)
    compute_timeseries!(sol.t, t₀)
end

function get_initial_conditions!(sol::SolutionODE{DT,TT}, q::Union{Vector{DT}, Vector{Double{DT}}}, k) where {DT,TT}
    get_data!(sol.q, q, 0, k)
end

function copy_solution!(sol::SolutionODE{DT,TT}, q::Union{Vector{DT}, Vector{Double{DT}}}, n, k) where {DT,TT}
    if mod(n, sol.nsave) == 0
        set_data!(sol.q, q, div(n, sol.nsave), k)
        sol.counter += 1
    end
end

function reset!(s::SolutionODE)
    reset!(s.q)
    compute_timeseries!(solution.t, solution.t[end])
    s.counter = 0
end


"Creates HDF5 file and initialises datasets for ODE solution object."
function create_hdf5(solution::SolutionODE{DT,TT,2}, file::AbstractString, ntime::Int=1) where {DT,TT}
    @assert ntime ≥ 1

    # create HDF5 file and save attributes and common parameters
    h5 = createHDF5(solution, file)

    # create dataset
    # ntime can be used to set the expected total number of timesteps
    # so that the size of the array does not need to be adapted dynamically.
    # Right now, it has to be set as dynamical size adaptation is not yet
    # working.
    q = d_create(h5, "q", datatype(DT), dataspace(solution.nd, solution.nt+1), "chunk", (solution.nd,1))

    # copy initial conditions
    q[1:solution.nd, 1] = solution.q.d[1:solution.nd, 1]

    return h5
end

"Creates HDF5 file and initialises datasets for ODE solution object."
function create_hdf5(solution::SolutionODE{DT,TT,3}, file::AbstractString, ntime::Int=1) where {DT,TT}
    @assert ntime ≥ 1

    # create HDF5 file and save attributes and common parameters
    h5 = createHDF5(solution, file)

    # create dataset
    # ntime can be used to set the expected total number of timesteps
    # so that the size of the array does not need to be adapted dynamically.
    # Right now, it has to be set as dynamical size adaptation is not yet
    # working.
    q = d_create(h5, "q", datatype(DT), dataspace(solution.nd, solution.nt+1, solution.ni), "chunk", (solution.nd,1,1))

    # copy initial conditions
    q[1:solution.nd, 1, 1:solution.ni] = solution.q.d[1:solution.nd, 1, 1:solution.ni]

    return h5
end

"Append solution to HDF5 file."
function CommonFunctions.write_to_hdf5(solution::SolutionODE{DT,TT,2}, h5::HDF5.HDF5File, offset=0) where {DT,TT}
    # set convenience variables and compute ranges
    d  = solution.nd
    n  = solution.nt
    j1 = offset+2
    j2 = offset+1+n

    # # extend dataset if necessary
    # if size(x, 2) < j2
    #     set_dims!(x, (d, j2))
    # end

    # copy data from solution to HDF5 dataset
    h5["q"][1:d, j1:j2] = solution.q.d[1:d, 2:n+1]

    return nothing
end

"Append solution to HDF5 file."
function CommonFunctions.write_to_hdf5(solution::SolutionODE{DT,TT,3}, h5::HDF5.HDF5File, offset=0) where {DT,TT}
    # set convenience variables and compute ranges
    d  = solution.nd
    n  = solution.nt
    i  = solution.ni
    j1 = offset+2
    j2 = offset+1+n

    # # extend dataset if necessary
    # if size(x, 2) < j2
    #     set_dims!(x, (d, j2))
    # end

    # copy data from solution to HDF5 dataset
    h5["q"][1:d, j1:j2, 1:i] = solution.q.d[1:d, 2:n+1, 1:i]

    return nothing
end
