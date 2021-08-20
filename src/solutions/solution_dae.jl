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

### Constructors

```julia
SSolutionDAE(equation, Δt, ntimesteps; nsave=DEFAULT_NSAVE, nwrite=DEFAULT_NWRITE, filename=nothing)
SSolutionDAE(t::TimeSeries, q::DataSeries, λ::DataSeries, ntimesteps)
SSolutionDAE(file::String)
PSolutionDAE(equation, Δt, ntimesteps; nsave=DEFAULT_NSAVE, nwrite=DEFAULT_NWRITE, filename=nothing)
PSolutionDAE(t::TimeSeries, q::PDataSeries, λ::PDataSeries, ntimesteps)
PSolutionDAE(file::String)
```

The constructors `SSolutionDAE` create a `SolutionDAE` with internal data structures
for serial simulations (i.e., standard arrays), while the constructors `PSolutionDAE`
create a `SolutionDAE` with internal data structures for parallel simulations (i.e.,
shared arrays).

The usual way to initialise a `Solution` is by passing an equation, which for
`SolutionDAE` has to be an [`DAE`](@ref), a time step `Δt` and the number of
time steps `ntimesteps`. The optional parameters `nsave` and `nwrite` determine
the intervals for storing the solution and writing to file, i.e., if
`nsave > 1` only every `nsave`'th solution is actually stored, and every
`nwrite`'th time step the solution is stored to disk.

The other constructors, either passing a `TimeSeries` and two `DataSeries` or a
filename are used to read data from previous simulations.

"""
abstract type SolutionDAE{dType, tType, N} <: DeterministicSolution{dType, tType, N} end

# Create SolutionPDAEs with serial and parallel data structures.
for (TSolution, TDataSeries, Tdocstring) in
    ((:SSolutionDAE, :DataSeries, "Serial Solution of a differential algebraic equation."),
    #  (:PSolutionDAE, :PDataSeries, "Parallel Solution of a differential algebraic equation.")
    )
    @eval begin
        $Tdocstring
        mutable struct $TSolution{dType, tType, N} <: SolutionDAE{dType, tType, N}
            nd::Int
            nm::Int
            nt::Int
            ni::Int
            t::TimeSeries{tType}
            q::$TDataSeries{dType,N}
            λ::$TDataSeries{dType,N}
            ntime::Int
            nsave::Int
            nwrite::Int
            counter::Vector{Int}
            woffset::Int
            periodicity::dType
            h5::HDF5.File

            function $TSolution{dType, tType, N}(nd, nm, nt, ni, t, q, λ, ntime, nsave, nwrite, periodicity=zero(q[begin])) where {dType <: Union{Number,AbstractArray}, tType <: Real, N}
                new(nd, nm, nt, ni, t, q, λ, ntime, nsave, nwrite, zeros(Int, ni), 0, periodicity)
            end
        end

        function $TSolution(equation::DAE{DT,TT1,AT}, Δt::TT2, ntimesteps::Int;
                            nsave::Int=DEFAULT_NSAVE, nwrite::Int=DEFAULT_NWRITE, filename=nothing) where {DT,TT1,TT2,AT}
            @assert nsave > 0
            @assert ntimesteps == 0 || ntimesteps ≥ nsave
            @assert nwrite == 0 || nwrite ≥ nsave
            @assert mod(ntimesteps, nsave) == 0

            if nwrite > 0
                @assert mod(nwrite, nsave) == 0
                @assert mod(ntimesteps, nwrite) == 0
            end

            TT = promote_type(TT1,TT2)
            N  = nsamples(equation) > 1 ? 2 : 1
            nd = ndims(equation)
            nm = equation.m
            ni = nsamples(equation)
            nt = div(ntimesteps, nsave)
            nt = (nwrite == 0 ? nt : div(nwrite, nsave))
            nw = (nwrite == 0 ? ntimesteps : nwrite)

            @assert nd > 0
            @assert nm > 0
            @assert ni > 0
            @assert nt ≥ 0
            @assert nw ≥ 0

            t = TimeSeries{TT}(nt, Δt, nsave)
            q = $TDataSeries(equation.q₀, nt, ni)
            λ = $TDataSeries(equation.λ₀, nt, ni)
            s = $TSolution{AT,TT,N}(nd, nm, nt, ni, t, q, λ, ntimesteps, nsave, nw, periodicity(equation))
            set_initial_conditions!(s, equation)

            if !isnothing(filename)
                isfile(filename) ? @warn("Overwriting existing HDF5 file.") : nothing
                create_hdf5!(s, filename)
            end

            return s
        end

        function $TSolution(t::TimeSeries{TT}, q::$TDataSeries{DT,N}, λ::$TDataSeries{DT,N}, ntimesteps::Int) where {DT,TT,N}
            @assert ndims(q) >= ndims(λ)
            @assert ntime(q) == ntime(λ)
            @assert nsamples(q) == nsamples(λ)

            # extract parameters
            nd = length(q[begin])
            nm = length(λ[begin])
            ni = nsamples(q)
            nt = t.n
            ns = div(ntimesteps, nt)

            @assert mod(ntimesteps, nt) == 0

            # create solution
            $TSolution{DT,TT,N}(nd, nm, nt, ni, t, q, λ, ntimesteps, ns, 0)
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
            λ = fromarray($TDataSeries, read(h5["λ"]), nsamples)

            # need to close the file
            close(h5)

            # create solution
            $TSolution(t, q, λ, ntimesteps)
        end
    end
end


Base.:(==)(sol1::SolutionDAE{DT1,TT1,N1}, sol2::SolutionDAE{DT2,TT2,N2}) where {DT1,TT1,N1,DT2,TT2,N2} = (
                                DT1 == DT2
                             && TT1 == TT2
                             && N1  == N2
                             && sol1.nd == sol2.nd
                             && sol1.nm == sol2.nm
                             && sol1.nt == sol2.nt
                             && sol1.ni == sol2.ni
                             && sol1.t  == sol2.t
                             && sol1.q  == sol2.q
                             && sol1.λ  == sol2.λ
                             && sol1.ntime == sol2.ntime
                             && sol1.nsave == sol2.nsave
                             && sol1.nwrite == sol2.nwrite
                             && sol1.counter == sol2.counter
                             && sol1.woffset == sol2.woffset
                             && sol1.periodicity == sol2.periodicity)

@inline hdf5(sol::SolutionDAE)  = sol.h5
@inline GeometricBase.counter(sol::SolutionDAE) = sol.counter
@inline GeometricBase.offset(sol::SolutionDAE) = sol.woffset
@inline GeometricBase.lastentry(sol::SolutionDAE) = sol.ni == 1 ? sol.counter[1] - 1 : sol.counter .- 1
@inline GeometricBase.nsamples(sol::SolutionDAE) = sol.ni
@inline GeometricBase.nsave(sol::SolutionDAE) = sol.nsave
@inline GeometricBase.ntime(sol::SolutionDAE) = sol.ntime
@inline GeometricBase.timesteps(sol::SolutionDAE)  = sol.t
@inline GeometricBase.eachtimestep(sol::SolutionDAE) = 1:sol.nt*sol.nsave
@inline GeometricBase.periodicity(sol::SolutionDAE) = sol.periodicity


function set_initial_conditions!(sol::SolutionDAE, equ::DAE)
    set_initial_conditions!(sol, equ.t₀, equ.q₀, equ.λ₀)
end

function set_initial_conditions!(sol::SolutionDAE{DT}, t₀::Real, q₀::DT, λ₀::DT) where {DT}
    set_data!(sol.q, q₀, 0)
    set_data!(sol.λ, λ₀, 0)
    compute_timeseries!(sol.t, t₀)
    sol.counter .= 1
end

function set_initial_conditions!(sol::SolutionDAE{DT}, t₀::Real, q₀::AbstractVector{DT}, λ₀::AbstractVector{DT}) where {DT}
    for i in eachindex(q₀,λ₀)
        set_data!(sol.q, q₀[i], 0, i)
        set_data!(sol.λ, λ₀[i], 0, i)
    end
    compute_timeseries!(sol.t, t₀)
    sol.counter .= 1
end

function get_initial_conditions!(sol::SolutionDAE{AT,TT}, asol::AtomicSolutionDAE{DT,TT,AT}, k, n=1) where {DT, TT, AT <: AbstractArray{DT}}
    get_solution!(sol, asol.q, asol.λ, n-1, k)
    asol.t  = sol.t[n-1]
    asol.q̃ .= 0
end

function get_initial_conditions!(sol::SolutionDAE{DT}, q::DT, λ::DT, k, n=1) where {DT}
    get_data!(sol.q, q, n-1, k)
    get_data!(sol.λ, λ, n-1, k)
end

function get_initial_conditions(sol::SolutionDAE, k, n=1)
    get_solution(sol, n-1, k)
end

function set_solution!(sol::SolutionDAE, t, q, λ, n, k=1)
    set_solution!(sol, q, λ, n, k)
end

function set_solution!(sol::SolutionDAE{AT,TT}, asol::AtomicSolutionDAE{DT,TT,AT}, n, k=1) where {DT, TT, AT <: AbstractArray{DT}}
    set_solution!(sol, asol.t, asol.q, asol.λ, n, k)
end

function set_solution!(sol::SolutionDAE{AT}, q::AT, λ::AT, n, k=1) where {AT}
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

function get_solution!(sol::SolutionDAE{AT}, q::AT, λ::AT, n, k=1) where {AT}
    get_data!(sol.q, q, n, k)
    get_data!(sol.λ, λ, n, k)
end

function get_solution(sol::SolutionDAE{AT,TT,1}, n, k=1) where {AT,TT}
    @assert k == 1
    (sol.t[n], sol.q[n], sol.λ[n])
end

function get_solution(sol::SolutionDAE{AT,TT,2}, n, k=1) where {AT,TT}
    (sol.t[n], sol.q[n,k], sol.λ[n,k])
end

function GeometricBase.reset!(sol::SolutionDAE)
    reset!(sol.q)
    reset!(sol.λ)
    compute_timeseries!(sol.t, sol.t[end])
    sol.counter .= 1
    sol.woffset += sol.nt
end

function Base.close(solution::SolutionDAE)
    close(solution.h5)
end
