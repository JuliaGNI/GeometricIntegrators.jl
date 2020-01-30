"""
`SolutionPDAE`: Solution of a partitioned differential algebraic equation

Contains all fields necessary to store the solution of an PDAE.

### Fields

* `nd`: dimension of the dynamical variable ``q``
* `nm`: dimension of the constraint submanifold
* `nt`: number of time steps to store
* `ni`: number of initial conditions
* `t`:  time steps
* `q`:  solution `q[nd, nt+1, ni]` with `q[:,0,:]` the initial conditions
* `p`:  solution `p[nd, nt+1, ni]` with `p[:,0,:]` the initial conditions
* `λ`:  Lagrange multiplier `λ[nd, nt+1, ni]`
* `ntime`: number of time steps to compute
* `nsave`: store every nsave'th time step (default: 1)
* `nwrite`: save data to disk after every nwrite'th time step (default: ntime)
* `counter`: counter for copied solution entries
* `woffset`: counter for file offset
* `h5`: HDF5 file for storage
"""
abstract type SolutionPDAE{dType, tType, N} <: DeterministicSolution{dType, tType, N} end

# Create SolutionPDAEs with serial and parallel data structures.
for (TSolution, TDataSeries, Tdocstring) in
    ((:SSolutionPDAE, :SDataSeries, "Serial Solution of a partitioned differential algebraic equation."),
     (:PSolutionPDAE, :PDataSeries, "Parallel Solution of a partitioned differential algebraic equation."))
    @eval begin
        $Tdocstring
        mutable struct $TSolution{dType, tType, N} <: SolutionPDAE{dType, tType, N}
            nd::Int
            nm::Int
            nt::Int
            ni::Int
            t::TimeSeries{tType}
            q::$TDataSeries{dType,N}
            p::$TDataSeries{dType,N}
            λ::$TDataSeries{dType,N}
            ntime::Int
            nsave::Int
            nwrite::Int
            counter::Vector{Int}
            woffset::Int
            h5::HDF5File

            function $TSolution{dType, tType, N}(nd, nm, nt, ni, t, q, p, λ, ntime, nsave, nwrite) where {dType <: Number, tType <: Real, N}
                new(nd, nm, nt, ni, t, q, p, λ, ntime, nsave, nwrite, zeros(Int, ni), 0)
            end
        end

        function $TSolution(equation::Union{IODE{DT,TT},VODE{DT,TT},PDAE{DT,TT},IDAE{DT,TT},VDAE{DT,TT}}, Δt::TT, ntime::Int, nsave::Int=DEFAULT_NSAVE, nwrite::Int=DEFAULT_NWRITE; filename=nothing) where {DT,TT}
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
            q = $TDataSeries(DT, nd, nt, ni)
            p = $TDataSeries(DT, nd, nt, ni)
            λ = $TDataSeries(DT, nm, nt, ni)
            s = $TSolution{DT,TT,N}(nd, nm, nt, ni, t, q, p, λ, ntime, nsave, nw)
            set_initial_conditions!(s, equation)

            if !isnothing(filename)
                isfile(filename) ? @warn("Overwriting existing HDF5 file.") : nothing
                create_hdf5(s, filename)
            end

            return s
        end

        function $TSolution(t::TimeSeries{TT}, q::$TDataSeries{DT,N}, p::$TDataSeries{DT,N}, λ::$TDataSeries{DT,N}, ntime::Int) where {DT,TT,N}
            @assert q.nd == p.nd >= λ.nd
            @assert q.nt == p.nt == λ.nt
            @assert q.ni == p.ni == λ.ni

            # extract parameters
            nd = q.nd
            nm = λ.nd
            ni = q.ni
            nt = t.n
            ns = div(ntime, nt)

            @assert mod(ntime, nt) == 0

            # create solution
            $TSolution{DT,TT,N}(nd, nm, nt, ni, t, q, p, λ, ntime, ns, 0)
        end

        function $TSolution(file::String)
            # open HDF5 file
            get_config(:verbosity) > 1 ? @info("Reading HDF5 file ", file) : nothing
            h5 = h5open(file, "r")

            # read attributes
            ntime = read(attrs(h5)["ntime"])
            nsave = read(attrs(h5)["nsave"])

            # reading data arrays
            t = TimeSeries(read(h5["t"]), nsave)
            q = $TDataSeries(read(h5["q"]))
            p = $TDataSeries(read(h5["p"]))
            λ = $TDataSeries(read(h5["λ"]))

            # need to close the file
            close(h5)

            # create solution
            $TSolution(t, q, p, λ, ntime)
        end
    end
end

hdf5(sol::SolutionPDAE)  = sol.h5
time(sol::SolutionPDAE)  = sol.t.t
ntime(sol::SolutionPDAE) = sol.ntime
nsave(sol::SolutionPDAE) = sol.nsave
offset(sol::SolutionPDAE) = sol.woffset


function set_initial_conditions!(sol::SolutionPDAE, equ::Union{IODE,VODE,PDAE,IDAE,VDAE})
    set_initial_conditions!(sol, equ.t₀, equ.q₀, equ.p₀, equ.λ₀)
end

function set_initial_conditions!(sol::SolutionPDAE{DT,TT}, t₀::TT, q₀::Union{Array{DT}, Array{TwicePrecision{DT}}}, p₀::Union{Array{DT}, Array{TwicePrecision{DT}}}, λ₀::Union{Array{DT}, Array{TwicePrecision{DT}}}) where {DT,TT}
    set_data!(sol.q, q₀, 0)
    set_data!(sol.p, p₀, 0)
    set_data!(sol.λ, λ₀, 0)
    compute_timeseries!(sol.t, t₀)
    sol.counter .= 1
end

function get_initial_conditions!(sol::SolutionPDAE{DT,TT}, asol::AtomicSolutionPDAE{DT,TT}, k, n=1) where {DT,TT}
    get_solution!(sol, asol.q, asol.p, asol.λ, n-1, k)
    asol.t  = sol.t[n-1]
    asol.q̃ .= 0
    asol.p̃ .= 0
end

function get_initial_conditions!(sol::SolutionPDAE{DT,TT}, q::SolutionVector{DT}, p::SolutionVector{DT}, λ::SolutionVector{DT}, k, n=1) where {DT,TT}
    get_solution!(sol, q, p, λ, n-1, k)
end

function get_initial_conditions!(sol::SolutionPDAE{DT,TT}, q::SolutionVector{DT}, p::SolutionVector{DT}, k, n=1) where {DT,TT}
    get_solution!(sol, q, p, n-1, k)
end

function get_initial_conditions(sol::SolutionPDAE, k, n=1)
    get_solution(sol, n-1, k)
end

function get_solution!(sol::SolutionPDAE{DT,TT}, q::SolutionVector{DT}, p::SolutionVector{DT}, λ::SolutionVector{DT}, n, k) where {DT,TT}
    q .= sol.q[:, n, k]
    p .= sol.p[:, n, k]
    λ .= sol.λ[:, n, k]
end

function get_solution!(sol::SolutionPDAE{DT,TT}, q::SolutionVector{DT}, p::SolutionVector{DT}, n, k) where {DT,TT}
    q .= sol.q[:, n, k]
    p .= sol.p[:, n, k]
end

function get_solution(sol::SolutionPDAE, n, k)
    (sol.t[n], sol.q[:, n, k], sol.p[:, n, k], sol.λ[:, n, k])
end

function set_solution!(sol::SolutionPDAE, t, q, p, n, k)
    set_solution!(sol, q, p, n, k)
end

function set_solution!(sol::SolutionPDAE, t, q, p, λ, n, k)
    set_solution!(sol, q, p, λ, n, k)
end

function set_solution!(sol::SolutionPDAE{DT,TT}, asol::AtomicSolutionPDAE{DT,TT}, n, k) where {DT,TT}
    set_solution!(sol, asol.t, asol.q, asol.p, asol.λ, n, k)
end

function set_solution!(sol::SolutionPDAE{DT,TT}, q::SolutionVector{DT}, p::SolutionVector{DT}, λ::SolutionVector{DT}, n, k) where {DT,TT}
    @assert n <= sol.ntime
    @assert k <= sol.ni
    if mod(n, sol.nsave) == 0
        if sol.counter[k] > sol.nt
            @error("Solution overflow. Call write_to_hdf5() and reset!() before continuing the simulation.")
        end
        set_data!(sol.q, q, sol.counter[k], k)
        set_data!(sol.p, p, sol.counter[k], k)
        set_data!(sol.λ, λ, sol.counter[k], k)
        sol.counter[k] += 1
    end
end

function set_solution!(sol::SolutionPDAE{DT,TT}, q::SolutionVector{DT}, p::SolutionVector{DT}, n, k) where {DT,TT}
    @assert n <= sol.ntime
    @assert k <= sol.ni
    if mod(n, sol.nsave) == 0
        if sol.counter[k] > sol.nt
            @error("Solution overflow. Call write_to_hdf5() and reset!() before continuing the simulation.")
        end
        set_data!(sol.q, q, sol.counter[k], k)
        set_data!(sol.p, p, sol.counter[k], k)
        sol.counter[k] += 1
    end
end

function CommonFunctions.reset!(sol::SolutionPDAE)
    reset!(sol.q)
    reset!(sol.p)
    reset!(sol.λ)
    compute_timeseries!(sol.t, sol.t[end])
    sol.counter .= 1
    sol.woffset += sol.nt
end

function Base.close(solution::SolutionPDAE)
    # close(solution.t)
    # close(solution.q)
    # close(solution.p)
    # close(solution.λ)
    close(solution.h5)
end


"Creates HDF5 file and initialises datasets for PDAE solution object with single initial condition."
function create_hdf5(solution::SolutionPDAE{DT,TT,2}, file::AbstractString) where {DT,TT}
    # create HDF5 file and save attributes and common parameters
    solution.h5 = createHDF5(solution, file)

    # save attributes
    save_attributes(solution)

    # create datasets
    nt = div(solution.ntime, solution.nsave)
    t = d_create(solution.h5, "t", datatype(TT), dataspace((nt+1,)), "chunk", (1,))
    q = d_create(solution.h5, "q", datatype(DT), dataspace(solution.nd, nt+1), "chunk", (solution.nd,1))
    p = d_create(solution.h5, "p", datatype(DT), dataspace(solution.nd, nt+1), "chunk", (solution.nd,1))
    λ = d_create(solution.h5, "λ", datatype(DT), dataspace(solution.nm, nt+1), "chunk", (solution.nm,1))

    # copy initial conditions
    t[1] = solution.t[0]
    q[:, 1] = solution.q[:, 0]
    p[:, 1] = solution.p[:, 0]
    λ[:, 1] = solution.λ[:, 0]

    return solution.h5
end

"Creates HDF5 file and initialises datasets for PDAE solution object with multiple initial conditions."
function create_hdf5(solution::SolutionPDAE{DT,TT,3}, file::AbstractString) where {DT,TT}
    # create HDF5 file and save attributes and common parameters
    solution.h5 = createHDF5(solution, file)

    # save attributes
    save_attributes(solution, solution.h5)

    # create datasets
    nt = div(solution.ntime, solution.nsave)
    t = d_create(solution.h5, "t", datatype(TT), dataspace((nt+1,)), "chunk", (1,))
    q = d_create(solution.h5, "q", datatype(DT), dataspace(solution.nd, nt+1, solution.ni), "chunk", (solution.nd,1,1))
    p = d_create(solution.h5, "p", datatype(DT), dataspace(solution.nd, nt+1, solution.ni), "chunk", (solution.nd,1,1))
    λ = d_create(solution.h5, "λ", datatype(DT), dataspace(solution.nm, nt+1, solution.ni), "chunk", (solution.nm,1,1))

    # copy initial conditions
    t[1] = solution.t[0]
    q[:, 1, :] = solution.q[:, 0, :]
    p[:, 1, :] = solution.p[:, 0, :]
    λ[:, 1, :] = solution.λ[:, 0, :]

    return solution.h5
end

"Append solution to HDF5 file."
function CommonFunctions.write_to_hdf5(solution::SolutionPDAE{DT,TT,2}, h5::HDF5File, offset=0) where {DT,TT}
    # set convenience variables and compute ranges
    n  = solution.nt
    j1 = offset+2
    j2 = offset+1+n

    # copy data from solution to HDF5 dataset
    h5["q"][:, j1:j2] = solution.q[:, 1:n]
    h5["p"][:, j1:j2] = solution.p[:, 1:n]
    h5["λ"][:, j1:j2] = solution.λ[:, 1:n]

    return nothing
end

"Append solution to HDF5 file."
function CommonFunctions.write_to_hdf5(solution::SolutionPDAE{DT,TT,3}, h5::HDF5File, offset=0) where {DT,TT}
    # set convenience variables and compute ranges
    n  = solution.nt
    j1 = offset+2
    j2 = offset+1+n

    # copy data from solution to HDF5 dataset
    h5["q"][:, j1:j2, :] = solution.q[:, 1:n, :]
    h5["p"][:, j1:j2, :] = solution.p[:, 1:n, :]
    h5["λ"][:, j1:j2, :] = solution.λ[:, 1:n, :]

    return nothing
end
