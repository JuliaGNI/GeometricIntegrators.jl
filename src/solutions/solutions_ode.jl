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
* `nsave`: store every nsave'th time step (default: 1)
* `nwrite`: save data to disk after every nwrite'th time step (default: ntime)
* `counter`: counter for copied solution entries
* `woffset`: counter for file offset
* `h5`: HDF5 file for storage
"""
mutable struct SolutionODE{dType, tType, N} <: DeterministicSolution{dType, tType, N}
    nd::Int
    nt::Int
    ni::Int
    t::TimeSeries{tType}
    q::SDataSeries{dType,N}
    ntime::Int
    nsave::Int
    nwrite::Int
    counter::Vector{Int}
    woffset::Int
    h5::HDF5File

    function SolutionODE{dType, tType, N}(nd, nt, ni, t, q, ntime, nsave, nwrite) where {dType <: Number, tType <: Real, N}
        new(nd, nt, ni, t, q, ntime, nsave, nwrite, zeros(Int, ni), 0)
    end
end

function SolutionODE(equation::Union{ODE{DT,TT,FT},SODE{DT,TT,FT}}, Δt::TT, ntime::Int, nsave::Int=DEFAULT_NSAVE, nwrite::Int=DEFAULT_NWRITE; filename=nothing) where {DT,TT,FT}
    @assert nsave > 0
    @assert ntime == 0 || ntime ≥ nsave
    @assert nwrite == 0 || nwrite ≥ nsave
    @assert mod(ntime, nsave) == 0

    if nwrite > 0
        @assert mod(nwrite, nsave) == 0
        @assert mod(ntime, nwrite) == 0
    end

    N  = (equation.n > 1 ? 3 : 2)
    nd = equation.d
    ni = equation.n
    nt = div(ntime, nsave)
    nt = (nwrite == 0 ? nt : div(nwrite, nsave))
    nw = (nwrite == 0 ? ntime : nwrite)

    @assert nd > 0
    @assert ni > 0
    @assert nt > 0
    @assert nw > 0

    t = TimeSeries{TT}(nt, Δt, nsave)
    q = SDataSeries(DT, nd, nt, ni)
    s = SolutionODE{DT,TT,N}(nd, nt, ni, t, q, ntime, nsave, nw)
    set_initial_conditions!(s, equation)

    if !isnothing(filename)
        isfile(filename) ? @warn("Overwriting existing HDF5 file.") : nothing
        create_hdf5(s, filename)
    end

    return s
end

function SolutionODE(t::TimeSeries{TT}, q::SDataSeries{DT,N}, ntime::Int) where {DT,TT,N}
    # extract parameters
    nd = q.nd
    ni = q.ni
    nt = t.n
    ns = div(ntime, nt)

    @assert mod(ntime, nt) == 0

    # create solution
    SolutionODE{DT,TT,N}(nd, nt, ni, t, q, ntime, ns, 0)
end

function SolutionODE(file::String)
    # open HDF5 file
    @info("Reading HDF5 file ", file)
    h5 = h5open(file, "r")

    # read attributes
    ntime = read(attrs(h5)["ntime"])
    nsave = read(attrs(h5)["nsave"])

    # reading data arrays
    t = TimeSeries(read(h5["t"]), nsave)
    q = SDataSeries(read(h5["q"]))

    # need to close the file
    close(h5)

    # create solution
    SolutionODE(t, q, ntime)
end

Base.:(==)(sol1::SolutionODE, sol2::SolutionODE) = (
                                sol1.nd == sol2.nd
                             && sol1.nt == sol2.nt
                             && sol1.ni == sol2.ni
                             && sol1.t  == sol2.t
                             && sol1.q  == sol2.q
                             && sol1.ntime == sol2.ntime
                             && sol1.nsave == sol2.nsave
                             && sol1.counter == sol2.counter)

hdf5(sol::SolutionODE)  = sol.h5
time(sol::SolutionODE)  = sol.t.t
ntime(sol::SolutionODE) = sol.ntime
nsave(sol::SolutionODE) = sol.nsave
offset(sol::SolutionODE) = sol.woffset


function set_initial_conditions!(sol::SolutionODE{DT,TT}, equ::Union{ODE{DT,TT},SODE{DT,TT}}) where {DT,TT}
    set_initial_conditions!(sol, equ.t₀, equ.q₀)
end

function set_initial_conditions!(sol::SolutionODE{DT,TT}, t₀::TT, q₀::Union{Array{DT}, Array{TwicePrecision{DT}}}) where {DT,TT}
    set_data!(sol.q, q₀, 0)
    compute_timeseries!(sol.t, t₀)
    sol.counter .= 1
end

function get_initial_conditions!(sol::SolutionODE{DT,TT}, q::Union{Vector{DT}, Vector{TwicePrecision{DT}}}, k) where {DT,TT}
    get_data!(sol.q, q, 0, k)
end

function get_initial_conditions(sol::SolutionODE, k, n=1)
    (sol.t[n], sol.q[:, n-1, k])
end

function CommonFunctions.get_solution!(sol::SolutionODE{DT,TT}, q::Union{Vector{DT}, Vector{TwicePrecision{DT}}}, n, k) where {DT,TT}
    q .= sol.q[:, n, k]
end

function CommonFunctions.get_solution(sol::SolutionODE, n, k)
    (sol.t[n], sol.q[:, n, k])
end

function CommonFunctions.set_solution!(sol, t, q, n, k)
    set_solution!(sol, q, n, k)
end

function CommonFunctions.set_solution!(sol::SolutionODE{DT,TT}, q::Union{Vector{DT}, Vector{TwicePrecision{DT}}}, n, k) where {DT,TT}
    @assert n <= sol.ntime
    @assert k <= sol.ni
    if mod(n, sol.nsave) == 0
        if sol.counter[k] > sol.nt
            @error("Solution overflow. Call write_to_hdf5() and reset!() before continuing the simulation.")
        end
        set_data!(sol.q, q, sol.counter[k], k)
        sol.counter[k] += 1
    end
end

function CommonFunctions.reset!(sol::SolutionODE)
    reset!(sol.q)
    compute_timeseries!(sol.t, sol.t[end])
    sol.counter .= 1
    sol.woffset += sol.nt
end

function Base.close(solution::SolutionODE)
    # close(solution.t)
    # close(solution.q)
    close(solution.h5)
end


"Creates HDF5 file and initialises datasets for ODE solution object."
function create_hdf5(solution::SolutionODE{DT,TT,2}, file::AbstractString) where {DT,TT}
    # create HDF5 file and save attributes and common parameters
    solution.h5 = createHDF5(solution, file)

    # save attributes
    save_attributes(solution)

    # create dataset
    nt = div(solution.ntime, solution.nsave)

    # nt can be used to set the expected total number of timesteps
    # so that the size of the array does not need to be adapted dynamically.
    # Right now, the size has to be set to nt as dynamical size adaptation is
    # not yet working.
    t = d_create(solution.h5, "t", datatype(TT), dataspace((nt+1,)), "chunk", (1,))
    q = d_create(solution.h5, "q", datatype(DT), dataspace(solution.nd, solution.nt+1), "chunk", (solution.nd,1))

    # copy initial conditions
    t[1] = solution.t[0]
    q[:, 1] = solution.q[:, 0]

    return solution.h5
end

"Creates HDF5 file and initialises datasets for ODE solution object."
function create_hdf5(solution::SolutionODE{DT,TT,3}, file::AbstractString) where {DT,TT}
    # create HDF5 file and save attributes and common parameters
    solution.h5 = createHDF5(solution, file)

    # create dataset
    nt = div(solution.ntime, solution.nsave)

    # nt can be used to set the expected total number of timesteps
    # so that the size of the array does not need to be adapted dynamically.
    # Right now, the size has to be set to nt as dynamical size adaptation is
    # not yet working.
    t = d_create(solution.h5, "t", datatype(TT), dataspace((nt+1,)), "chunk", (1,))
    q = d_create(solution.h5, "q", datatype(DT), dataspace(solution.nd, nt+1, solution.ni), "chunk", (solution.nd,1,1))

    # copy initial conditions
    t[1] = solution.t[0]
    q[:, 1, :] = solution.q[:, 0, :]

    return solution.h5
end

"Append solution to HDF5 file."
function CommonFunctions.write_to_hdf5(solution::SolutionODE{DT,TT,2}, h5::HDF5File, offset=0) where {DT,TT}
    # set convenience variables and compute ranges
    n  = solution.nt
    j1 = offset+2
    j2 = offset+1+n

    # TODO # extend dataset if necessary
    # if size(x, 2) < j2
    #     set_dims!(x, (d, j2))
    # end

    # copy data from solution to HDF5 dataset
    h5["t"][j1:j2] = solution.t[1:n]
    h5["q"][:, j1:j2] = solution.q[:, 1:n]

    return nothing
end

"Append solution to HDF5 file."
function CommonFunctions.write_to_hdf5(solution::SolutionODE{DT,TT,3}, h5::HDF5File, offset=0) where {DT,TT}
    # set convenience variables and compute ranges
    n  = solution.nt
    j1 = offset+2
    j2 = offset+1+n

    # TODO # extend dataset if necessary
    # if size(x, 2) < j2
    #     set_dims!(x, (d, j2))
    # end

    # copy data from solution to HDF5 dataset
    h5["t"][j1:j2] = solution.t[1:n]
    h5["q"][:, j1:j2, :] = solution.q[:, 1:n, :]

    return nothing
end
