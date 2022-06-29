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
* `offset`: counter for file offset
* `counter`: counter for copied solution entries

### Constructors

```julia
SolutionODE(equation, Δt, ntimesteps; nsave=DEFAULT_NSAVE, nwrite=DEFAULT_NWRITE, filename=nothing)
SolutionODE(t::TimeSeries, q::DataSeries, ntimesteps)
SolutionODE(h5io::SolutionHDF5)
SolutionODE(file::String)
```

The usual way to initialise a `Solution` is by passing an equation, which for
`SolutionODE` has to be an [`ODE`](@ref) or [`SODE`](@ref), a time step `Δt`
and the number of time steps `ntimesteps`. The optional parameters `nsave` and
`nwrite` determine the intervals for storing the solution and writing to file,
i.e., if `nsave > 1` only every `nsave`'th solution is actually stored, and
every `nwrite`'th time step the solution is stored to disk.

The other constructors, either passing a `TimeSeries` and a `DataSeries` or a
filename are used to read data from previous simulations.

"""
mutable struct SolutionODE{dType,tType,N} <: DeterministicSolution{dType,tType,N}
    nd::Int
    nt::Int
    ni::Int
    t::TimeSeries{tType}
    q::DataSeries{dType,N}
    ntime::Int
    nsave::Int
    nwrite::Int
    offset::Int
    counter::Vector{Int}
    periodicity::dType

    function SolutionODE{dType,tType,N}(nd, nt, ni, t, q, ntime, nsave, nwrite, periodicity = NullPeriodicity()) where {dType<:Union{Number,AbstractArray},tType<:Real,N}
        new(nd, nt, ni, t, q, ntime, nsave, nwrite, 0, zeros(Int, ni), _periodicity(q[0], periodicity))
    end
end

function SolutionODE(equation::Union{ODEProblem{DT,TT1,AT},SODEProblem{DT,TT1,AT}}, Δt::TT2, ntimesteps::Int;
    nsave::Int = DEFAULT_NSAVE, nwrite::Int = DEFAULT_NWRITE) where {DT,TT1,TT2,AT}
    @assert nsave > 0
    @assert ntimesteps == 0 || ntimesteps ≥ nsave
    @assert nwrite == 0 || nwrite ≥ nsave
    @assert mod(ntimesteps, nsave) == 0

    if nwrite > 0
        @assert mod(nwrite, nsave) == 0
        @assert mod(ntimesteps, nwrite) == 0
    end

    TT = promote_type(TT1, TT2)
    N = nsamples(equation) > 1 ? 2 : 1
    nd = length(vec(equation.ics.q))
    ni = nsamples(equation)
    nt = div(ntimesteps, nsave)
    nt = (nwrite == 0 ? nt : div(nwrite, nsave))
    nw = (nwrite == 0 ? ntimesteps : nwrite)

    @assert nd > 0
    @assert ni > 0
    @assert nt ≥ 0
    @assert nw ≥ 0

    t = TimeSeries{TT}(nt, Δt, nsave)
    q = DataSeries(equation.ics.q, nt, ni)
    s = SolutionODE{AT,TT,N}(nd, nt, ni, t, q, ntimesteps, nsave, nw, periodicity(equation))
    set_initial_conditions!(s, equation)

    return s
end

function SolutionODE(t::TimeSeries{TT}, q::DataSeries{DT,N}, ntimesteps::Int) where {DT,TT,N}
    # extract parameters
    nd = length(q[begin])
    ni = nsamples(q)
    nt = t.n
    ns = div(ntimesteps, nt)

    @assert mod(ntimesteps, nt) == 0

    # create solution
    SolutionODE{DT,TT,N}(nd, nt, ni, t, q, ntimesteps, ns, 0)
end

function SolutionODE(h5io::SolutionHDF5)
    # get hdf5 file descriptor
    h5 = hdf5(h5io)

    # read attributes
    ntimesteps = read(attributes(h5)["ntime"])
    nsave = read(attributes(h5)["nsave"])
    nsamples = read(attributes(h5)["nsamples"])

    # reading data arrays
    t = TimeSeries(read(h5["t"]), nsave)
    q = fromarray(DataSeries, read(h5["q"]), nsamples)

    # create solution
    SolutionODE(t, q, ntimesteps)
end

function SolutionODE(file::String)
    load(file) do h5io
        SolutionODE(h5io::SolutionHDF5)
    end
end


Base.:(==)(sol1::SolutionODE{DT1,TT1,N1}, sol2::SolutionODE{DT2,TT2,N2}) where {DT1,TT1,N1,DT2,TT2,N2} = (
    DT1 == DT2
    && TT1 == TT2
    && N1 == N2
    && sol1.nd == sol2.nd
    && sol1.nt == sol2.nt
    && sol1.ni == sol2.ni
    && sol1.t == sol2.t
    && sol1.q == sol2.q
    && sol1.ntime == sol2.ntime
    && sol1.nsave == sol2.nsave
    && sol1.nwrite == sol2.nwrite
    && sol1.offset == sol2.offset
    && sol1.counter == sol2.counter
    && sol1.periodicity == sol2.periodicity)

@inline GeometricBase.counter(sol::SolutionODE) = sol.counter
@inline offset(sol::SolutionODE) = sol.offset
@inline GeometricBase.lastentry(sol::SolutionODE) = sol.ni == 1 ? sol.counter[1] - 1 : sol.counter .- 1
@inline GeometricBase.nsamples(sol::SolutionODE) = sol.ni
@inline GeometricBase.nsave(sol::SolutionODE) = sol.nsave
@inline GeometricBase.ntime(sol::SolutionODE) = sol.ntime
@inline GeometricBase.timesteps(sol::SolutionODE) = sol.t
@inline GeometricBase.eachtimestep(sol::SolutionODE) = 1:sol.nt*sol.nsave
@inline GeometricBase.periodicity(sol::SolutionODE) = sol.periodicity


function set_initial_conditions!(sol::SolutionODE, equ::AbstractProblemODE)
    set_initial_conditions!(sol, equ.tspan[begin], equ.ics.q)
end

function set_initial_conditions!(sol::SolutionODE{DT}, t₀::Real, q₀) where {DT}
    set_data!(sol.q, q₀, 0)
    compute_timeseries!(sol.t, t₀)
    sol.counter .= 1
end

function set_initial_conditions!(sol::SolutionODE{DT}, t₀::Real, q₀::AbstractVector{DT}) where {DT}
    for i in eachindex(q₀)
        set_data!(sol.q, q₀[i], 0, i)
    end
    compute_timeseries!(sol.t, t₀)
    sol.counter .= 1
end

function get_initial_conditions!(sol::SolutionODE{AT,TT}, asol::AtomicSolutionODE{DT,TT,AT}, k, n = 1) where {DT,TT,AT<:AbstractArray{DT}}
    get_solution!(sol, asol.q, n - 1, k)
    asol.t = sol.t[n-1]
    asol.q̃ .= 0
end

function get_initial_conditions!(sol::SolutionODE{DT}, q::DT, k, n = 1) where {DT}
    get_solution!(sol, q, n - 1, k)
end

function get_initial_conditions(sol::SolutionODE, k, n = 1)
    get_solution(sol, n - 1, k)
end

function set_solution!(sol::SolutionODE, t, q, n, k = 1)
    set_solution!(sol, q, n, k)
end

function set_solution!(sol::SolutionODE{AT,TT}, asol::AtomicSolutionODE{DT,TT,AT}, n, k = 1) where {DT,TT,AT<:AbstractArray{DT}}
    set_solution!(sol, asol.t, asol.q, n, k)
end

function set_solution!(sol::SolutionODE{AT}, q::AT, n, k = 1) where {AT}
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

function get_solution!(sol::SolutionODE{AT}, q::AT, n, k = 1) where {AT}
    get_data!(sol.q, q, n, k)
end

function get_solution(sol::SolutionODE{AT,TT,1}, n, k = 1) where {AT,TT}
    @assert k == 1
    (sol.t[n], sol.q[n])
end

function get_solution(sol::SolutionODE{AT,TT,2}, n, k = 1) where {AT,TT}
    (sol.t[n], sol.q[n, k])
end

function GeometricBase.reset!(sol::SolutionODE)
    reset!(sol.q)
    compute_timeseries!(sol.t, sol.t[end])
    sol.counter .= 1
    sol.offset += sol.nt
end

function Base.close(solution::SolutionODE)
    close(solution.h5)
end
