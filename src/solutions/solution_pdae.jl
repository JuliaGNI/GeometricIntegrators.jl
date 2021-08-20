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

### Constructors

```julia
SSolutionPDAE(equation, Δt, ntimesteps; nsave=DEFAULT_NSAVE, nwrite=DEFAULT_NWRITE, filename=nothing)
SSolutionPDAE(t::TimeSeries, q::DataSeries, p::DataSeries, λ::DataSeries, ntimesteps)
SSolutionPDAE(file::String)
PSolutionPDAE(equation, Δt, ntimesteps; nsave=DEFAULT_NSAVE, nwrite=DEFAULT_NWRITE, filename=nothing)
PSolutionPDAE(t::TimeSeries, q::PDataSeries, p::PDataSeries, λ::PDataSeries, ntimesteps)
PSolutionPDAE(file::String)
```

The constructors `SSolutionPDAE` create a `SolutionPDAE` with internal data structures
for serial simulations (i.e., standard arrays), while the constructors `PSolutionPDAE`
create a `SolutionPDAE` with internal data structures for parallel simulations (i.e.,
shared arrays).

The usual way to initialise a `Solution` is by passing an equation, which for
`SolutionPDAE` has to be an [`PDAE`](@ref), [`HDAE`](@ref), [`IDAE`](@ref) or
[`LDAE`](@ref), a time step `Δt` and the number of time steps `ntimesteps`. The
optional parameters `nsave` and `nwrite` determine the intervals for storing
the solution and writing to file, i.e., if `nsave > 1` only every `nsave`'th
solution is actually stored, and every `nwrite`'th time step the solution is
stored to disk.

The other constructors, either passing a `TimeSeries` and three `DataSeries` or a
filename are used to read data from previous simulations.

"""
abstract type SolutionPDAE{dType, tType, N} <: DeterministicSolution{dType, tType, N} end

# Create SolutionPDAEs with serial and parallel data structures.
for (TSolution, TDataSeries, Tdocstring) in
    ((:SSolutionPDAE, :DataSeries, "Serial Solution of a partitioned differential algebraic equation."),
    #  (:PSolutionPDAE, :PDataSeries, "Parallel Solution of a partitioned differential algebraic equation.")
    )
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
            periodicity::dType
            h5::HDF5.File

            function $TSolution{dType, tType, N}(nd, nm, nt, ni, t, q, p, λ, ntime, nsave, nwrite, periodicity=zero(q[begin])) where {dType <: Union{Number,AbstractArray}, tType <: Real, N}
                new(nd, nm, nt, ni, t, q, p, λ, ntime, nsave, nwrite, zeros(Int, ni), 0, periodicity)
            end
        end

        function $TSolution(equation::Union{PDAE{DT,TT1,AT}, HDAE{DT,TT1,AT}, IDAE{DT,TT1,AT}, LDAE{DT,TT1,AT}, IODE{DT,TT1,AT}, LODE{DT,TT1,AT}}, Δt::TT2, ntimesteps::Int;
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
            p = $TDataSeries(equation.p₀, nt, ni)
            λ = $TDataSeries(equation.λ₀, nt, ni)
            s = $TSolution{AT,TT,N}(nd, nm, nt, ni, t, q, p, λ, ntimesteps, nsave, nw, periodicity(equation))
            set_initial_conditions!(s, equation)

            if !isnothing(filename)
                isfile(filename) ? @warn("Overwriting existing HDF5 file.") : nothing
                create_hdf5!(s, filename)
            end

            return s
        end

        function $TSolution(t::TimeSeries{TT}, q::$TDataSeries{DT,N}, p::$TDataSeries{DT,N}, λ::$TDataSeries{DT,N}, ntimesteps::Int) where {DT,TT,N}
            @assert ndims(q) == ndims(p) >= ndims(λ)
            @assert ntime(q) == ntime(p) == ntime(λ)
            @assert nsamples(q) == nsamples(p) == nsamples(λ)

            # extract parameters
            nd = length(q[begin])
            nm = length(λ[begin])
            ni = nsamples(q)
            nt = t.n
            ns = div(ntimesteps, nt)

            @assert mod(ntimesteps, nt) == 0

            # create solution
            $TSolution{DT,TT,N}(nd, nm, nt, ni, t, q, p, λ, ntimesteps, ns, 0)
        end

        function $TSolution(file::String)
            # open HDF5 file
            get_config(:verbosity) > 1 ? @info("Reading HDF5 file ", file) : nothing
            h5 = h5open(file, "r")

            # read attributes
            ntime = read(attributes(h5)["ntime"])
            nsave = read(attributes(h5)["nsave"])
            nsamples = read(attributes(h5)["nsamples"])

            # reading data arrays
            t = TimeSeries(read(h5["t"]), nsave)
            q = fromarray($TDataSeries, read(h5["q"]), nsamples)
            p = fromarray($TDataSeries, read(h5["p"]), nsamples)
            λ = fromarray($TDataSeries, read(h5["λ"]), nsamples)

            # need to close the file
            close(h5)

            # create solution
            $TSolution(t, q, p, λ, ntime)
        end
    end
end


Base.:(==)(sol1::SolutionPDAE{DT1,TT1,N1}, sol2::SolutionPDAE{DT2,TT2,N2}) where {DT1,TT1,N1,DT2,TT2,N2} = (
                                DT1 == DT2
                             && TT1 == TT2
                             && N1  == N2
                             && sol1.nd == sol2.nd
                             && sol1.nm == sol2.nm
                             && sol1.nt == sol2.nt
                             && sol1.ni == sol2.ni
                             && sol1.t  == sol2.t
                             && sol1.q  == sol2.q
                             && sol1.p  == sol2.p
                             && sol1.λ  == sol2.λ
                             && sol1.ntime == sol2.ntime
                             && sol1.nsave == sol2.nsave
                             && sol1.nwrite == sol2.nwrite
                             && sol1.counter == sol2.counter
                             && sol1.woffset == sol2.woffset
                             && sol1.periodicity == sol2.periodicity)

@inline hdf5(sol::SolutionPDAE)  = sol.h5
@inline GeometricBase.counter(sol::SolutionPDAE) = sol.counter
@inline GeometricBase.offset(sol::SolutionPDAE) = sol.woffset
@inline GeometricBase.lastentry(sol::SolutionPDAE) = sol.ni == 1 ? sol.counter[1] - 1 : sol.counter .- 1
@inline GeometricBase.nsamples(sol::SolutionPDAE) = sol.ni
@inline GeometricBase.nsave(sol::SolutionPDAE) = sol.nsave
@inline GeometricBase.ntime(sol::SolutionPDAE) = sol.ntime
@inline GeometricBase.timesteps(sol::SolutionPDAE)  = sol.t
@inline GeometricBase.eachtimestep(sol::SolutionPDAE) = 1:sol.nt*sol.nsave
@inline GeometricBase.periodicity(sol::SolutionPDAE) = sol.periodicity


function set_initial_conditions!(sol::SolutionPDAE, equ::Union{IODE,LODE,PDAE,IDAE,LDAE})
    set_initial_conditions!(sol, equ.t₀, equ.q₀, equ.p₀, equ.λ₀)
end

function set_initial_conditions!(sol::SolutionPDAE{DT}, t₀::Real, q₀::DT, p₀::DT, λ₀::DT) where {DT}
    set_data!(sol.q, q₀, 0)
    set_data!(sol.p, p₀, 0)
    set_data!(sol.λ, λ₀, 0)
    compute_timeseries!(sol.t, t₀)
    sol.counter .= 1
end

function set_initial_conditions!(sol::SolutionPDAE{DT}, t₀::Real, q₀::AbstractVector{DT}, p₀::AbstractVector{DT}, λ₀::AbstractVector{DT}) where {DT}
    for i in eachindex(q₀,p₀,λ₀)
        set_data!(sol.q, q₀[i], 0, i)
        set_data!(sol.p, p₀[i], 0, i)
        set_data!(sol.λ, λ₀[i], 0, i)
    end
    compute_timeseries!(sol.t, t₀)
    sol.counter .= 1
end

function get_initial_conditions!(sol::SolutionPDAE{AT,TT}, asol::AtomicSolutionPDAE{DT,TT,AT}, k, n=1) where {DT,TT,AT}
    get_solution!(sol, asol.q, asol.p, asol.λ, n-1, k)
    asol.t  = sol.t[n-1]
    asol.q̃ .= 0
    asol.p̃ .= 0
end

function get_initial_conditions!(sol::SolutionPDAE{DT}, q::DT, p::DT, λ::DT, k, n=1) where {DT}
    get_data!(sol.q, q, n-1, k)
    get_data!(sol.p, p, n-1, k)
    get_data!(sol.λ, λ, n-1, k)
end

function get_initial_conditions!(sol::SolutionPDAE{DT}, q::DT, p::DT, k, n=1) where {DT}
    get_data!(sol.q, q, n-1, k)
    get_data!(sol.p, p, n-1, k)
end

function get_initial_conditions(sol::SolutionPDAE, k, n=1)
    get_solution(sol, n-1, k)
end

function set_solution!(sol::SolutionPDAE, t::Real, q::Vector, p::Vector, n::Int, k::Int=1)
    set_solution!(sol, q, p, n, k)
end

function set_solution!(sol::SolutionPDAE, t::Real, q::Vector, p::Vector, λ::Vector, n::Int, k::Int=1)
    set_solution!(sol, q, p, λ, n, k)
end

function set_solution!(sol::SolutionPDAE{AT,TT}, asol::AtomicSolutionPDAE{DT,TT,AT}, n, k=1) where {DT,TT,AT}
    set_solution!(sol, asol.t, asol.q, asol.p, asol.λ, n, k)
end

function set_solution!(sol::SolutionPDAE{AT}, q::AT, p::AT, λ::AT, n, k=1) where {AT}
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

function set_solution!(sol::SolutionPDAE{AT}, q::AT, p::AT, n, k=1) where {AT}
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

function get_solution!(sol::SolutionPDAE{AT}, q::AT, p::AT, λ::AT, n, k=1) where {AT}
    get_data!(sol.q, q, n, k)
    get_data!(sol.p, p, n, k)
    get_data!(sol.λ, λ, n, k)
end

function get_solution!(sol::SolutionPDAE{AT}, q::AT, p::AT, n, k=1) where {AT}
    get_data!(sol.q, q, n, k)
    get_data!(sol.p, p, n, k)
end

function get_solution(sol::SolutionPDAE{AT,TT,1}, n, k=1) where {AT,TT}
    @assert k == 1
    (sol.t[n], sol.q[n], sol.p[n], sol.λ[n])
end

function get_solution(sol::SolutionPDAE{AT,TT,2}, n, k=1) where {AT,TT}
    (sol.t[n], sol.q[n,k], sol.p[n,k], sol.λ[n,k])
end

function GeometricBase.reset!(sol::SolutionPDAE)
    reset!(sol.q)
    reset!(sol.p)
    reset!(sol.λ)
    compute_timeseries!(sol.t, sol.t[end])
    sol.counter .= 1
    sol.woffset += sol.nt
end

function Base.close(solution::SolutionPDAE)
    close(solution.h5)
end
