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
* `offset`: counter for file offset
* `counter`: counter for copied solution entries

### Constructors

```julia
SolutionPODE(equation, Δt, ntimesteps; nsave=DEFAULT_NSAVE, nwrite=DEFAULT_NWRITE, filename=nothing)
SolutionPODE(t::TimeSeries, q::DataSeries, p::DataSeries, ntimesteps)
SolutionPODE(h5io::SolutionHDF5)
SolutionPODE(file::String)
```

The constructors `SSolutionPODE` create a `SolutionPODE` with internal data structures
for serial simulations (i.e., standard arrays), while the constructors `PSolutionPODE`
create a `SolutionPODE` with internal data structures for parallel simulations (i.e.,
shared arrays).

The usual way to initialise a `Solution` is by passing an equation, which for
`SolutionPODE` has to be an [`PODE`](@ref), [`HODE`](@ref), [`IODE`](@ref) or
[`LODE`](@ref), a time step `Δt` and the number of time steps `ntimesteps`.
The optional parameters `nsave` and `nwrite` determine the intervals for
storing the solution and writing to file, i.e., if `nsave > 1` only every
`nsave`'th solution is actually stored, and every `nwrite`'th time step the
solution is stored to disk.

The other constructors, either passing a `TimeSeries` and two `DataSeries` or a
filename are used to read data from previous simulations.

"""
mutable struct SolutionPODE{dType,tType,N} <: DeterministicSolution{dType,tType,N}
    nd::Int
    nt::Int
    ni::Int
    t::TimeSeries{tType}
    q::DataSeries{dType,N}
    p::DataSeries{dType,N}
    ntime::Int
    nsave::Int
    nwrite::Int
    offset::Int
    counter::Vector{Int}
    periodicity::dType

    function SolutionPODE{dType,tType,N}(nd, nt, ni, t, q, p, ntime, nsave, nwrite, periodicity = zero(q[begin])) where {dType<:Union{Number,AbstractArray},tType<:Real,N}
        new(nd, nt, ni, t, q, p, ntime, nsave, nwrite, 0, zeros(Int, ni), periodicity)
    end
end

function SolutionPODE(equation::Union{PODE{DT,TT1,AT},HODE{DT,TT1,AT},IODE{DT,TT1,AT},LODE{DT,TT1,AT}}, Δt::TT2, ntimesteps::Int;
    nsave::Int = DEFAULT_NSAVE, nwrite::Int = DEFAULT_NWRITE, filename = nothing) where {DT,TT1,TT2,AT}
    @assert nsave > 0
    @assert ntimesteps == 0 || ntimesteps ≥ nsave
    @assert nwrite == 0 || nwrite ≥ nsave
    @assert mod(ntimesteps, nsave) == 0

    if nwrite > 0
        @assert mod(nwrite, nsave) == 0
        @assert mod(ntimesteps, nwrite) == 0
    end

    TT = promote_type(TT1, TT2)
    N = (nsamples(equation) > 1 ? 2 : 1)
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
    q = DataSeries(equation.q₀, nt, ni)
    p = DataSeries(equation.p₀, nt, ni)
    s = SolutionPODE{AT,TT,N}(nd, nt, ni, t, q, p, ntimesteps, nsave, nw, periodicity(equation))
    set_initial_conditions!(s, equation)

    if !isnothing(filename)
        isfile(filename) ? @warn("Overwriting existing HDF5 file.") : nothing
        create_hdf5!(s, filename)
    end

    return s
end

function SolutionPODE(t::TimeSeries{TT}, q::DataSeries{AT,N}, p::DataSeries{AT,N}, ntimesteps::Int) where {AT,TT,N}
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
    SolutionPODE{AT,TT,N}(nd, nt, ni, t, q, p, ntimesteps, ns, 0)
end

function SolutionPODE(h5io::SolutionHDF5)
    # get hdf5 file descriptor
    h5 = hdf5(h5io)

    # read attributes
    ntimesteps = read(attributes(h5)["ntime"])
    nsave = read(attributes(h5)["nsave"])
    nsamples = read(attributes(h5)["nsamples"])

    # reading data arrays
    t = TimeSeries(read(h5["t"]), nsave)
    q = fromarray(DataSeries, read(h5["q"]), nsamples)
    p = fromarray(DataSeries, read(h5["p"]), nsamples)

    # create solution
    SolutionPODE(t, q, p, ntimesteps)
end

function SolutionPODE(file::String)
    load(file) do h5io
        SolutionPODE(h5io::SolutionHDF5)
    end
end


Base.:(==)(sol1::SolutionPODE, sol2::SolutionPODE) = (
    sol1.nd == sol2.nd
    && sol1.nt == sol2.nt
    && sol1.ni == sol2.ni
    && sol1.t == sol2.t
    && sol1.q == sol2.q
    && sol1.p == sol2.p
    && sol1.ntime == sol2.ntime
    && sol1.nsave == sol2.nsave
    && sol1.nwrite == sol2.nwrite
    && sol1.offset == sol2.offset
    && sol1.counter == sol2.counter
    && sol1.periodicity == sol2.periodicity)

@inline GeometricBase.counter(sol::SolutionPODE) = sol.counter
@inline GeometricBase.offset(sol::SolutionPODE) = sol.offset
@inline GeometricBase.lastentry(sol::SolutionPODE) = sol.ni == 1 ? sol.counter[1] - 1 : sol.counter .- 1
@inline GeometricBase.nsamples(sol::SolutionPODE) = sol.ni
@inline GeometricBase.nsave(sol::SolutionPODE) = sol.nsave
@inline GeometricBase.ntime(sol::SolutionPODE) = sol.ntime
@inline GeometricBase.timesteps(sol::SolutionPODE) = sol.t
@inline GeometricBase.eachtimestep(sol::SolutionPODE) = 1:sol.nt*sol.nsave
@inline GeometricBase.periodicity(sol::SolutionPODE) = sol.periodicity


function set_initial_conditions!(sol::SolutionPODE, equ::Union{PODE,HODE,IODE,LODE})
    set_initial_conditions!(sol, equ.t₀, equ.q₀, equ.p₀)
end

function set_initial_conditions!(sol::SolutionPODE{AT}, t₀::Real, q₀::AT, p₀::AT) where {AT}
    set_data!(sol.q, q₀, 0)
    set_data!(sol.p, p₀, 0)
    compute_timeseries!(sol.t, t₀)
    sol.counter .= 1
end

function set_initial_conditions!(sol::SolutionPODE{AT}, t₀::Real, q₀::AbstractVector{AT}, p₀::AbstractVector{AT}) where {AT}
    for i in eachindex(q₀, p₀)
        set_data!(sol.q, q₀[i], 0, i)
        set_data!(sol.p, p₀[i], 0, i)
    end
    compute_timeseries!(sol.t, t₀)
    sol.counter .= 1
end

function get_initial_conditions!(sol::SolutionPODE{AT,TT}, asol::AtomicSolutionPODE{DT,TT,AT}, k, n = 1) where {DT,TT,AT<:AbstractArray{DT}}
    get_solution!(sol, asol.q, asol.p, n - 1, k)
    asol.t = sol.t[n-1]
    asol.q̃ .= 0
    asol.p̃ .= 0
end

function get_initial_conditions!(sol::SolutionPODE{AT}, q::AT, p::AT, k, n = 1) where {AT}
    get_solution!(sol, q, p, n - 1, k)
end

function get_initial_conditions(sol::SolutionPODE, k, n = 1)
    get_solution(sol, n - 1, k)
end

function get_solution!(sol::SolutionPODE{AT}, q::AT, p::AT, n, k = 1) where {AT}
    get_data!(sol.q, q, n, k)
    get_data!(sol.p, p, n, k)
end

function get_solution(sol::SolutionPODE{AT,TT,1}, n, k = 1) where {AT,TT}
    @assert k == 1
    (sol.t[n], sol.q[n], sol.p[n])
end

function get_solution(sol::SolutionPODE{AT,TT,2}, n, k = 1) where {AT,TT}
    (sol.t[n], sol.q[n, k], sol.p[n, k])
end

function set_solution!(sol::SolutionPODE, t, q, p, n, k)
    set_solution!(sol, q, p, n, k)
end

function set_solution!(sol::SolutionPODE{AT,TT}, asol::AtomicSolutionPODE{DT,TT,AT}, n, k = 1) where {DT,TT,AT<:AbstractArray{DT}}
    set_solution!(sol, asol.t, asol.q, asol.p, n, k)
end

function set_solution!(sol::SolutionPODE{AT,TT}, q::AT, p::AT, n, k = 1) where {AT,TT}
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

function GeometricBase.reset!(sol::SolutionPODE)
    reset!(sol.q)
    reset!(sol.p)
    compute_timeseries!(sol.t, sol.t[end])
    sol.counter .= 1
    sol.offset += sol.nt
end

function Base.close(solution::SolutionPODE)
    close(solution.h5)
end
