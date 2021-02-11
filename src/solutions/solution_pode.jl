"""
`SolutionPODE`: Solution of a partitioned ordinary differential equation

Contains all fields necessary to store the solution of an PODE.

### Fields

* `nd`: dimension of the dynamical variable ``q``
* `nt`: number of time steps to store
* `ni`: number of initial conditions
* `t`:  time steps
* `q`:  solution `q[nd, nt+1, ni]` with `q[:,0,:]` the initial conditions
* `p`:  solution `p[nd, nt+1, ni]` with `p[:,0,:]` the initial conditions
* `ntime`: number of time steps to compute
* `nsave`: store every nsave'th time step (default: 1)
* `nwrite`: save data to disk after every nwrite'th time step (default: ntime)
* `counter`: counter for copied solution entries
* `woffset`: counter for file offset
* `h5`: HDF5 file for storage

### Constructors

```julia
SSolutionPODE(equation, Δt, ntimesteps; nsave=DEFAULT_NSAVE, nwrite=DEFAULT_NWRITE, filename=nothing)
SSolutionPODE(t::TimeSeries, q::SDataSeries, p::SDataSeries, ntimesteps)
SSolutionPODE(file::String)
PSolutionPODE(equation, Δt, ntimesteps; nsave=DEFAULT_NSAVE, nwrite=DEFAULT_NWRITE, filename=nothing)
PSolutionPODE(t::TimeSeries, q::PDataSeries, p::PDataSeries, ntimesteps)
PSolutionPODE(file::String)
```

The constructors `SSolutionPODE` create a `SolutionPODE` with internal data structures
for serial simulations (i.e., standard arrays), while the constructors `PSolutionPODE`
create a `SolutionPODE` with internal data structures for parallel simulations (i.e.,
shared arrays).

The usual way to initialise a `Solution` is by passing an equation, which for
`SolutionPODE` has to be an [`PODE`](@ref), [`HODE`](@ref), [`IODE`](@ref) or
[`VODE`](@ref), a time step `Δt` and the number of time steps `ntimesteps`.
The optional parameters `nsave` and `nwrite` determine the intervals for
storing the solution and writing to file, i.e., if `nsave > 1` only every
`nsave`'th solution is actually stored, and every `nwrite`'th time step the
solution is stored to disk.

The other constructors, either passing a `TimeSeries` and two `DataSeries` or a
filename are used to read data from previous simulations.

"""
abstract type SolutionPODE{dType, tType, N} <: DeterministicSolution{dType, tType, N} end

# Create SolutionPODEs with serial and parallel data structures.
for (TSolution, TDataSeries, Tdocstring) in
    ((:SSolutionPODE, :SDataSeries, "Serial Solution of a partitioned ordinary differential equation."),
     (:PSolutionPODE, :PDataSeries, "Parallel Solution of a partitioned ordinary differential equation."))
    @eval begin
        $Tdocstring
        mutable struct $TSolution{dType, tType, N} <: SolutionPODE{dType, tType, N}
            nd::Int
            nt::Int
            ni::Int
            t::TimeSeries{tType}
            q::$TDataSeries{dType,N}
            p::$TDataSeries{dType,N}
            ntime::Int
            nsave::Int
            nwrite::Int
            counter::Vector{Int}
            woffset::Int
            periodicity::dType
            h5::HDF5.File

            function $TSolution{dType, tType, N}(nd, nt, ni, t, q, p, ntime, nsave, nwrite, periodicity=zero(q[begin])) where {dType <: Union{Number,AbstractArray}, tType <: Real, N}
                new(nd, nt, ni, t, q, p, ntime, nsave, nwrite, zeros(Int, ni), 0, periodicity)
            end
        end

        function $TSolution(equation::Union{PODE{DT,TT,AT}, HODE{DT,TT,AT}, IODE{DT,TT,AT}, VODE{DT,TT,AT}}, Δt::TT, ntimesteps::Int;
                            nsave::Int=DEFAULT_NSAVE, nwrite::Int=DEFAULT_NWRITE, filename=nothing) where {DT,TT,AT}
            @assert nsave > 0
            @assert ntimesteps == 0 || ntimesteps ≥ nsave
            @assert nwrite == 0 || nwrite ≥ nsave
            @assert mod(ntimesteps, nsave) == 0

            if nwrite > 0
                @assert mod(nwrite, nsave) == 0
                @assert mod(ntimesteps, nwrite) == 0
            end

            N  = (nsamples(equation) > 1 ? 2 : 1)
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
            p = $TDataSeries(equation.p₀, nt, ni)
            s = $TSolution{AT,TT,N}(nd, nt, ni, t, q, p, ntimesteps, nsave, nw, periodicity(equation))
            set_initial_conditions!(s, equation)

            if !isnothing(filename)
                isfile(filename) ? @warn("Overwriting existing HDF5 file.") : nothing
                create_hdf5!(s, filename)
            end

            return s
        end

        function $TSolution(t::TimeSeries{TT}, q::$TDataSeries{AT,N}, p::$TDataSeries{AT,N}, ntimesteps::Int) where {AT,TT,N}
            @assert ndims(q) == ndims(p)
            @assert ntime(q) == ntime(p)
            @assert nsamples(q) == nsamples(p)

            # extract parameters
            nd = ndims(q)
            ni = nsamples(q)
            nt = t.n
            ns = div(ntimesteps, nt)

            @assert mod(ntimesteps, nt) == 0

            # create solution
            $TSolution{AT,TT,N}(nd, nt, ni, t, q, p, ntimesteps, ns, 0)
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
            p = fromarray($TDataSeries, read(h5["p"]), nsamples)

            # need to close the file
            close(h5)

            # create solution
            $TSolution(t, q, p, ntimesteps)
        end
    end
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
                             && sol1.counter == sol2.counter
                             && sol1.periodicity == sol2.periodicity)

@inline hdf5(sol::SolutionPODE)  = sol.h5
@inline timesteps(sol::SolutionPODE)  = sol.t
@inline nsave(sol::SolutionPODE) = sol.nsave
@inline counter(sol::SolutionPODE) = sol.counter
@inline offset(sol::SolutionPODE) = sol.woffset
@inline lastentry(sol::SolutionPODE) = sol.ni == 1 ? sol.counter[1] - 1 : sol.counter .- 1
@inline Common.ntime(sol::SolutionPODE) = sol.ntime
@inline Common.periodicity(sol::SolutionPODE) = sol.periodicity


function set_initial_conditions!(sol::SolutionPODE, equ::Union{PODE,HODE,IODE,VODE})
    set_initial_conditions!(sol, equ.t₀, equ.q₀, equ.p₀)
end

function set_initial_conditions!(sol::SolutionPODE{AT,TT}, t₀::TT, q₀::AT, p₀::AT) where {AT,TT}
    set_data!(sol.q, q₀, 0)
    set_data!(sol.p, p₀, 0)
    compute_timeseries!(sol.t, t₀)
    sol.counter .= 1
end

function set_initial_conditions!(sol::SolutionPODE{AT,TT}, t₀::TT, q₀::AbstractVector{AT}, p₀::AbstractVector{AT}) where {AT,TT}
    for i in eachindex(q₀,p₀)
        set_data!(sol.q, q₀[i], 0, i)
        set_data!(sol.p, p₀[i], 0, i)
    end
    compute_timeseries!(sol.t, t₀)
    sol.counter .= 1
end

function get_initial_conditions!(sol::SolutionPODE{AT,TT}, asol::AtomicSolutionPODE{DT,TT,AT}, k, n=1) where {DT, TT, AT <: AbstractArray{DT}}
    get_solution!(sol, asol.q, asol.p, n-1, k)
    asol.t  = sol.t[n-1]
    asol.q̃ .= 0
    asol.p̃ .= 0
end

function get_initial_conditions!(sol::SolutionPODE{AT,TT}, q::AT, p::AT, k, n=1) where {AT,TT}
    get_solution!(sol, q, p, n-1, k)
end

function get_initial_conditions(sol::SolutionPODE, k, n=1)
    get_solution(sol, n-1, k)
end

function get_solution!(sol::SolutionPODE{AT,TT}, q::AT, p::AT, n, k=1) where {AT,TT}
    get_data!(sol.q, q, n, k)
    get_data!(sol.p, p, n, k)
end

function get_solution(sol::SolutionPODE{AT,TT,1}, n, k=1) where {AT,TT}
    @assert k == 1
    (sol.t[n], sol.q[n], sol.p[n])
end

function get_solution(sol::SolutionPODE{AT,TT,2}, n, k=1) where {AT,TT}
    (sol.t[n], sol.q[n,k], sol.p[n,k])
end

function set_solution!(sol::SolutionPODE, t, q, p, n, k)
    set_solution!(sol, q, p, n, k)
end

function set_solution!(sol::SolutionPODE{AT,TT}, asol::AtomicSolutionPODE{DT,TT,AT}, n, k=1) where {DT, TT, AT <: AbstractArray{DT}}
    set_solution!(sol, asol.t, asol.q, asol.p, n, k)
end

function set_solution!(sol::SolutionPODE{AT,TT}, q::AT, p::AT, n, k=1) where {AT,TT}
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

function Common.reset!(sol::SolutionPODE)
    reset!(sol.q)
    reset!(sol.p)
    compute_timeseries!(sol.t, sol.t[end])
    sol.counter .= 1
    sol.woffset += sol.nt
end

function Base.close(solution::SolutionPODE)
    close(solution.h5)
end
