"""
`SolutionDAE`: Solution of a differential algebraic equation

Contains all fields necessary to store the solution of an DAE.

### Fields

* `nd`: dimension of the dynamical variable ``q``
* `nm`: dimension of the constraint submanifold
* `nt`: number of time steps to store
* `ni`: number of initial conditions
* `t`:  time steps
* `q`:  solution `q[nd, nt+1, ni]` with `q[:,0,:]` the initial conditions
* `λ`:  Lagrange multiplier `λ[nd, nt+1, ni]`
* `ntime`: number of time steps to compute
* `nsave`: store every nsave'th time step (default: 1)
* `nwrite`: save data to disk after every nwrite'th time step (default: ntime)
* `counter`: counter for copied solution entries
* `woffset`: counter for file offset
* `h5`: HDF5 file for storage
"""
mutable struct SolutionDAE{dType, tType, N} <: DeterministicSolution{dType, tType, N}
    nd::Int
    nm::Int
    nt::Int
    ni::Int
    t::TimeSeries{tType}
    q::SDataSeries{dType,N}
    λ::SDataSeries{dType,N}
    ntime::Int
    nsave::Int
    nwrite::Int
    counter::Vector{Int}
    woffset::Int
    h5::HDF5File

    function SolutionDAE{dType, tType, N}(nd, nm, nt, ni, t, q, λ, ntime, nsave, nwrite) where {dType <: Number, tType <: Real, N}
        new(nd, nm, nt, ni, t, q, λ, ntime, nsave, nwrite, zeros(Int, ni), 0)
    end
end

function SolutionDAE(equation::DAE{DT,TT}, Δt::TT, ntime::Int, nsave::Int=DEFAULT_NSAVE, nwrite::Int=DEFAULT_NWRITE; filename=nothing) where {DT,TT}
    @assert nsave > 0
    @assert ntime == 0 || ntime ≥ nsave
    @assert nwrite == 0 || nwrite ≥ nsave
    @assert mod(ntime, nsave) == 0

    if nwrite > 0
        @assert mod(nwrite, nsave) == 0
        @assert mod(ntime, nwrite) == 0
    end

    N  = equation.n > 1 ? 3 : 2
    nd = equation.d
    nm = equation.m
    ni = equation.n
    nt = div(ntime, nsave)
    nt = (nwrite == 0 ? nt : div(nwrite, nsave))
    nw = (nwrite == 0 ? ntime : nwrite)

    @assert nd > 0
    @assert nm > 0
    @assert ni > 0
    @assert nw > 0

    t = TimeSeries{TT}(nt, Δt, nsave)
    q = SDataSeries(DT, nd, nt, ni)
    λ = SDataSeries(DT, nm, nt, ni)
    s = SolutionDAE{DT,TT,N}(nd, nm, nt, ni, t, q, λ, ntime, nsave, nw)
    set_initial_conditions!(s, equation)

    if !isnothing(filename)
        isfile(filename) ? @warn("Overwriting existing HDF5 file.") : nothing
        create_hdf5(s, filename)
    end

    return s
end

function SolutionDAE(t::TimeSeries{TT}, q::SDataSeries{DT,N}, λ::SDataSeries{DT,N}, ntime::Int) where {DT,TT,N}
    @assert q.nd >= λ.nd
    @assert q.nt == λ.nt
    @assert q.ni == λ.ni

    # extract parameters
    nd = q.nd
    nm = λ.nd
    ni = q.ni
    nt = t.n
    ns = div(ntime, nt)

    @assert mod(ntime, nt) == 0

    # create solution
    SolutionDAE{DT,TT,N}(nd, nm, nt, ni, t, q, λ, ntime, ns, 0)
end

function SolutionDAE(file::String)
    # open HDF5 file
    get_config(:verbosity) > 1 ? @info("Reading HDF5 file ", file) : nothing
    h5 = h5open(file, "r")

    # read attributes
    ntime = read(attrs(h5)["ntime"])
    nsave = read(attrs(h5)["nsave"])

    # reading data arrays
    t = TimeSeries(read(h5["t"]), nsave)
    q = SDataSeries(read(h5["q"]))
    λ = SDataSeries(read(h5["λ"]))

    # need to close the file
    close(h5)

    # create solution
    SolutionDAE(t, q, λ, ntime)
end


hdf5(sol::SolutionDAE)  = sol.h5
time(sol::SolutionDAE)  = sol.t.t
ntime(sol::SolutionDAE) = sol.ntime
nsave(sol::SolutionDAE) = sol.nsave
offset(sol::SolutionDAE) = sol.woffset


function set_initial_conditions!(sol::SolutionDAE{DT,TT}, equ::DAE{DT,TT}) where {DT,TT}
    set_initial_conditions!(sol, equ.t₀, equ.q₀, equ.λ₀)
end

function set_initial_conditions!(sol::SolutionDAE{DT,TT}, t₀::TT, q₀::Array{DT}, λ₀::Array{DT}) where {DT,TT}
    set_data!(sol.q, q₀, 0)
    set_data!(sol.λ, λ₀, 0)
    compute_timeseries!(sol.t, t₀)
    sol.counter .= 1
end

function get_initial_conditions!(sol::SolutionDAE{DT,TT}, asol::AtomisticSolutionDAE{DT,TT}, k, n=1) where {DT,TT}
    get_solution!(sol, asol.q, asol.λ, n-1, k)
    asol.t  = sol.t[n-1]
    asol.q̃ .= 0
end

function get_initial_conditions!(sol::SolutionDAE{DT,TT}, q::Vector{DT}, λ::Vector{DT}, k) where {DT,TT}
    get_data!(sol.q, q, 0, k)
    get_data!(sol.λ, λ, 0, k)
end

function CommonFunctions.set_solution!(sol::SolutionDAE, t, q, λ, n, k)
    set_solution!(sol, q, λ, n, k)
end

function CommonFunctions.set_solution!(sol::SolutionDAE{DT,TT}, asol::AtomisticSolutionDAE{DT,TT}, n, k) where {DT,TT}
    set_solution!(sol, asol.t, asol.q, asol.λ, n, k)
end

function CommonFunctions.set_solution!(sol::SolutionDAE{DT,TT}, q::Vector{DT}, λ::Vector{DT}, n, k) where {DT,TT}
    @assert n <= sol.ntime
    @assert k <= sol.ni
    if mod(n, sol.nsave) == 0
        if sol.counter[k] > sol.nt
            @error("Solution overflow. Call write_to_hdf5() and reset!() before continuing the simulation.")
        end
        set_data!(sol.q, q, sol.counter[k], k)
        set_data!(sol.λ, λ, sol.counter[k], k)
        sol.counter[k] += 1
    end
end

function CommonFunctions.reset!(sol::SolutionDAE)
    reset!(sol.q)
    reset!(sol.λ)
    compute_timeseries!(sol.t, sol.t[end])
    sol.counter .= 1
    sol.woffset += sol.nt
end

function Base.close(solution::SolutionDAE)
    # close(solution.t)
    # close(solution.q)
    # close(solution.λ)
    close(solution.h5)
end


"Creates HDF5 file and initialises datasets for DAE solution object."
function create_hdf5(solution::SolutionDAE{DT,TT,2}, file::AbstractString) where {DT,TT}
    # create HDF5 file and save attributes and common parameters
    solution.h5 = createHDF5(solution, file)

    # save attributes
    save_attributes(solution)

    # create datasets
    nt = div(solution.ntime, solution.nsave)
    t = d_create(solution.h5, "t", datatype(TT), dataspace((nt+1,)), "chunk", (1,))
    q = d_create(solution.h5, "q", datatype(DT), dataspace(solution.nd, nt+1), "chunk", (solution.nd,1))
    λ = d_create(solution.h5, "λ", datatype(DT), dataspace(solution.nm, nt+1), "chunk", (solution.nm,1))

    # copy initial conditions
    t[1] = solution.t[0]
    q[:, 1] = solution.q[:, 0]
    λ[:, 1] = solution.λ[:, 0]

    return solution.h5
end

"Creates HDF5 file and initialises datasets for DAE solution object."
function create_hdf5(solution::SolutionDAE{DT,TT,3}, file::AbstractString) where {DT,TT}
    # create HDF5 file and save attributes and common parameters
    solution.h5 = createHDF5(solution, file)

    # save attributes
    save_attributes(solution)

    # create datasets
    nt = div(solution.ntime, solution.nsave)
    t = d_create(solution.h5, "t", datatype(TT), dataspace((nt+1,)), "chunk", (1,))
    q = d_create(solution.h5, "q", datatype(DT), dataspace(solution.nd, nt+1, solution.ni), "chunk", (solution.nd,1,1))
    λ = d_create(solution.h5, "λ", datatype(DT), dataspace(solution.nm, nt+1, solution.ni), "chunk", (solution.nm,1,1))

    # copy initial conditions
    t[1] = solution.t[0]
    q[:, 1, :] = solution.q[:, 0, :]
    λ[:, 1, :] = solution.λ[:, 0, :]

    return solution.h5
end

"Append solution to HDF5 file."
function CommonFunctions.write_to_hdf5(solution::SolutionDAE{DT,TT,2}, h5::HDF5File, offset=0) where {DT,TT}
    # set convenience variables and compute ranges
    n  = solution.nt
    j1 = offset+2
    j2 = offset+1+n

    # copy data from solution to HDF5 dataset
    h5["q"][:, j1:j2] = solution.q[:, 1:n]
    h5["λ"][:, j1:j2] = solution.λ[:, 1:n]

    return nothing
end

"Append solution to HDF5 file."
function CommonFunctions.write_to_hdf5(solution::SolutionDAE{DT,TT,3}, h5::HDF5File, offset=0) where {DT,TT}
    # set convenience variables and compute ranges
    n  = solution.nt
    j1 = offset+2
    j2 = offset+1+n

    # copy data from solution to HDF5 dataset
    h5["q"][:, j1:j2, :] = solution.q[:, 1:n, :]
    h5["λ"][:, j1:j2, :] = solution.λ[:, 1:n, :]

    return nothing
end
