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
            periodicity::dType
            h5::HDF5.File

            function $TSolution{dType, tType, N}(nd, nt, ni, t, q, ntime, nsave, nwrite, periodicity=zero(q[begin])) where {dType <: Union{Number,AbstractArray}, tType <: Real, N}
                new(nd, nt, ni, t, q, ntime, nsave, nwrite, zeros(Int, ni), 0, periodicity)
            end
        end

        function $TSolution(equation::Union{ODE{DT,TT,AT},SODE{DT,TT,AT}}, Δt::TT, ntimesteps::Int;
                            nsave::Int=DEFAULT_NSAVE, nwrite::Int=DEFAULT_NWRITE, filename=nothing) where {DT,TT,AT}
            @assert nsave > 0
            @assert ntimesteps == 0 || ntimesteps ≥ nsave
            @assert nwrite == 0 || nwrite ≥ nsave
            @assert mod(ntimesteps, nsave) == 0

            if nwrite > 0
                @assert mod(nwrite, nsave) == 0
                @assert mod(ntimesteps, nwrite) == 0
            end

            N  = nsamples(equation) > 1 ? 2 : 1
            nd = ndims(equation)
            ni = nsamples(equation)
            nt = div(ntimesteps, nsave)
            nt = (nwrite == 0 ? nt : div(nwrite, nsave))
            nw = (nwrite == 0 ? ntimesteps : nwrite)

            @assert nd > 0
            @assert ni > 0
            @assert nt ≥ 0
            @assert nw ≥ 0

            t = TimeSeries{TT}(nt, Δt, nsave)
            q = $TDataSeries(equation.q₀, nt, ni)
            s = $TSolution{AT,TT,N}(nd, nt, ni, t, q, ntimesteps, nsave, nw, periodicity(equation))
            set_initial_conditions!(s, equation)

            if !isnothing(filename)
                isfile(filename) ? @warn("Overwriting existing HDF5 file.") : nothing
                create_hdf5!(s, filename)
            end

            return s
        end

        function $TSolution(t::TimeSeries{TT}, q::$TDataSeries{DT,N}, ntimesteps::Int) where {DT,TT,N}
            # extract parameters
            nd = length(q[begin])
            ni = nsamples(q)
            nt = t.n
            ns = div(ntimesteps, nt)

            @assert mod(ntimesteps, nt) == 0

            # create solution
            $TSolution{DT,TT,N}(nd, nt, ni, t, q, ntimesteps, ns, 0)
        end

        function $TSolution(file::String)
            # open HDF5 file
            get_config(:verbosity) > 1 ? @info("Reading HDF5 file ", file) : nothing
            h5 = h5open(file, "r")

            # read attributes
            ntimesteps = read(attributes(h5)["ntime"])
            nsave = read(attributes(h5)["nsave"])
            nsamples = read(attributes(h5)["nsamples"])

            # reading data arrays
            t = TimeSeries(read(h5["t"]), nsave)
            q = fromarray($TDataSeries, read(h5["q"]), nsamples)

            # need to close the file
            close(h5)

            # create solution
            $TSolution(t, q, ntimesteps)
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
                             && sol1.woffset == sol2.woffset
                             && sol1.periodicity == sol2.periodicity)

@inline hdf5(sol::SolutionODE)  = sol.h5
@inline timesteps(sol::SolutionODE)  = sol.t
@inline nsave(sol::SolutionODE) = sol.nsave
@inline counter(sol::SolutionODE) = sol.counter
@inline offset(sol::SolutionODE) = sol.woffset
@inline lastentry(sol::SolutionODE) = sol.ni == 1 ? sol.counter[1] - 1 : sol.counter .- 1
@inline Common.ntime(sol::SolutionODE) = sol.ntime
@inline Common.periodicity(sol::SolutionODE) = sol.periodicity


function set_initial_conditions!(sol::SolutionODE, equ::Union{ODE,SODE})
    set_initial_conditions!(sol, equ.t₀, equ.q₀)
end

function set_initial_conditions!(sol::SolutionODE{DT,TT}, t₀::TT, q₀::DT) where {DT,TT}
    set_data!(sol.q, q₀, 0)
    compute_timeseries!(sol.t, t₀)
    sol.counter .= 1
end

function set_initial_conditions!(sol::SolutionODE{DT,TT}, t₀::TT, q₀::AbstractVector{DT}) where {DT,TT}
    for i in eachindex(q₀)
        set_data!(sol.q, q₀[i], 0, i)
    end
    compute_timeseries!(sol.t, t₀)
    sol.counter .= 1
end

function get_initial_conditions!(sol::SolutionODE{AT,TT}, asol::AtomicSolutionODE{DT,TT,AT}, k, n=1) where {DT, TT, AT <: AbstractArray{DT}}
    get_solution!(sol, asol.q, n-1, k)
    asol.t  = sol.t[n-1]
    asol.q̃ .= 0
end

function get_initial_conditions!(sol::SolutionODE{DT}, q::DT, k, n=1) where {DT}
    get_solution!(sol, q, n-1, k)
end

function get_initial_conditions(sol::SolutionODE, k, n=1)
    get_solution(sol, n-1, k)
end

function set_solution!(sol::SolutionODE, t, q, n, k=1)
    set_solution!(sol, q, n, k)
end

function set_solution!(sol::SolutionODE{AT,TT}, asol::AtomicSolutionODE{DT,TT,AT}, n, k=1) where {DT, TT, AT <: AbstractArray{DT}}
    set_solution!(sol, asol.t, asol.q, n, k)
end

function set_solution!(sol::SolutionODE{AT}, q::AT, n, k=1) where {AT}
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

function get_solution!(sol::SolutionODE{AT}, q::AT, n, k=1) where {AT}
    get_data!(sol.q, q, n, k)
end

function get_solution(sol::SolutionODE{AT,TT,1}, n, k=1) where {AT,TT}
    @assert k == 1
    (sol.t[n], sol.q[n])
end

function get_solution(sol::SolutionODE{AT,TT,2}, n, k=1) where {AT,TT}
    (sol.t[n], sol.q[n,k])
end

function Common.reset!(sol::SolutionODE)
    reset!(sol.q)
    compute_timeseries!(sol.t, sol.t[end])
    sol.counter .= 1
    sol.woffset += sol.nt
end

function Base.close(solution::SolutionODE)
    close(solution.h5)
end
