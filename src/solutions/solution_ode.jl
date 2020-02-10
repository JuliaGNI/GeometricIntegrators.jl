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
abstract type SolutionODE{dType, tType, N} <: DeterministicSolution{dType, tType, N} end

# Create SolutionODEs with serial and parallel data structures.
for (TSolution, TDataSeries, Tdocstring) in
    ((:SSolutionODE, :SDataSeries, "Serial Solution of a ordinary differential equation."),
     (:PSolutionODE, :PDataSeries, "Parallel Solution of a ordinary differential equation."))
    @eval begin
        $Tdocstring
        mutable struct $TSolution{dType, tType, N} <: SolutionODE{dType, tType, N}
            nd::Int
            nt::Int
            ni::Int
            t::TimeSeries{tType}
            q::$TDataSeries{dType,N}
            ntime::Int
            nsave::Int
            nwrite::Int
            counter::Vector{Int}
            woffset::Int
            h5::HDF5File

            function $TSolution{dType, tType, N}(nd, nt, ni, t, q, ntime, nsave, nwrite) where {dType <: Number, tType <: Real, N}
                new(nd, nt, ni, t, q, ntime, nsave, nwrite, zeros(Int, ni), 0)
            end
        end

        function $TSolution(equation::Union{ODE{DT,TT,FT},SODE{DT,TT,FT}}, Δt::TT, ntime::Int; nsave::Int=DEFAULT_NSAVE, nwrite::Int=DEFAULT_NWRITE, filename=nothing) where {DT,TT,FT}
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
            @assert nt ≥ 0
            @assert nw ≥ 0

            t = TimeSeries{TT}(nt, Δt, nsave)
            q = $TDataSeries(DT, nd, nt, ni)
            s = $TSolution{DT,TT,N}(nd, nt, ni, t, q, ntime, nsave, nw)
            set_initial_conditions!(s, equation)

            if !isnothing(filename)
                isfile(filename) ? @warn("Overwriting existing HDF5 file.") : nothing
                create_hdf5!(s, filename)
            end

            return s
        end

        function $TSolution(t::TimeSeries{TT}, q::$TDataSeries{DT,N}, ntime::Int) where {DT,TT,N}
            # extract parameters
            nd = q.nd
            ni = q.ni
            nt = t.n
            ns = div(ntime, nt)

            @assert mod(ntime, nt) == 0

            # create solution
            $TSolution{DT,TT,N}(nd, nt, ni, t, q, ntime, ns, 0)
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

            # need to close the file
            close(h5)

            # create solution
            $TSolution(t, q, ntime)
        end
    end
end


Base.:(==)(sol1::SolutionODE{DT1,TT1,N1}, sol2::SolutionODE{DT2,TT2,N2}) where {DT1,TT1,N1,DT2,TT2,N2} = (
                                DT1 == DT2
                             && TT1 == TT2
                             && N1  == N2
                             && sol1.nd == sol2.nd
                             && sol1.nt == sol2.nt
                             && sol1.ni == sol2.ni
                             && sol1.t  == sol2.t
                             && sol1.q  == sol2.q
                             && sol1.ntime == sol2.ntime
                             && sol1.nsave == sol2.nsave
                             && sol1.nwrite == sol2.nwrite
                             && sol1.counter == sol2.counter
                             && sol1.woffset == sol2.woffset)

hdf5(sol::SolutionODE)  = sol.h5
timesteps(sol::SolutionODE)  = sol.t
ntime(sol::SolutionODE) = sol.ntime
nsave(sol::SolutionODE) = sol.nsave
offset(sol::SolutionODE) = sol.woffset


function set_initial_conditions!(sol::SolutionODE, equ::Union{ODE,SODE})
    set_initial_conditions!(sol, equ.t₀, equ.q₀)
end

function set_initial_conditions!(sol::SolutionODE{DT,TT}, t₀::TT, q₀::Union{Array{DT}, Array{TwicePrecision{DT}}}) where {DT,TT}
    set_data!(sol.q, q₀, 0)
    compute_timeseries!(sol.t, t₀)
    sol.counter .= 1
end

function get_initial_conditions!(sol::SolutionODE{DT,TT}, asol::AtomicSolutionODE{DT,TT}, k, n=1) where {DT,TT}
    get_solution!(sol, asol.q, n-1, k)
    asol.t  = sol.t[n-1]
    asol.q̃ .= 0
end

function get_initial_conditions!(sol::SolutionODE{DT}, q::SolutionVector{DT}, k, n=1) where {DT}
    get_solution!(sol, q, n-1, k)
end

function get_initial_conditions(sol::SolutionODE, k, n=1)
    get_solution(sol, n-1, k)
end

function get_solution!(sol::SolutionODE{DT}, q::SolutionVector{DT}, n, k=1) where {DT}
    for i in eachindex(q) q[i] = sol.q[i, n, k] end
end

function get_solution(sol::SolutionODE, n, k=1)
    (sol.t[n], sol.q[:, n, k])
end

function set_solution!(sol::SolutionODE, t, q, n, k=1)
    set_solution!(sol, q, n, k)
end

function set_solution!(sol::SolutionODE{DT,TT}, asol::AtomicSolutionODE{DT,TT}, n, k=1) where {DT,TT}
    set_solution!(sol, asol.t, asol.q, n, k)
end

function set_solution!(sol::SolutionODE{DT}, q::SolutionVector{DT}, n, k=1) where {DT}
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
    close(solution.h5)
end
