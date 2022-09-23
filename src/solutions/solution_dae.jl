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
* `offset`: counter for file offset
* `counter`: counter for copied solution entries

### Constructors

```julia
SolutionDAE(equation, Δt, ntimesteps; nsave=DEFAULT_NSAVE, nwrite=DEFAULT_NWRITE, filename=nothing)
SolutionDAE(t::TimeSeries, q::DataSeries, λ::DataSeries, ntimesteps)
SolutionDAE(h5io::SolutionHDF5)
SolutionDAE(file::String)
```

The usual way to initialise a `Solution` is by passing an equation, which for
`SolutionDAE` has to be an [`DAE`](@ref), a time step `Δt` and the number of
time steps `ntimesteps`. The optional parameters `nsave` and `nwrite` determine
the intervals for storing the solution and writing to file, i.e., if
`nsave > 1` only every `nsave`'th solution is actually stored, and every
`nwrite`'th time step the solution is stored to disk.

The other constructors, either passing a `TimeSeries` and two `DataSeries` or a
filename are used to read data from previous simulations.

"""
# mutable struct SolutionDAE{dType,tType,N} <: DeterministicSolution{dType,tType,N}
#     nd::Int
#     nm::Int
#     nt::Int
#     ni::Int
#     t::TimeSeries{tType}
#     q::DataSeries{dType,N}
#     λ::DataSeries{dType,N}
#     ntime::Int
#     nsave::Int
#     nwrite::Int
#     offset::Int
#     counter::Vector{Int}
#     periodicity::dType

#     function SolutionDAE{dType,tType,N}(nd, nm, nt, ni, t, q, λ, ntime, nsave, nwrite, periodicity = NullPeriodicity()) where {dType<:Union{Number,AbstractArray},tType<:Real,N}
#         new(nd, nm, nt, ni, t, q, λ, ntime, nsave, nwrite, 0, zeros(Int, ni), _periodicity(q[0], periodicity))
#     end
# end

# function SolutionDAE(problem::DAEProblem{DT,TT,AT};
#                      nsave::Int = DEFAULT_NSAVE, nwrite::Int = DEFAULT_NWRITE) where {DT,TT,AT}

#     ntimesteps = ntime(problem)
    
#     @assert nsave > 0
#     @assert ntimesteps == 0 || ntimesteps ≥ nsave
#     @assert nwrite == 0 || nwrite ≥ nsave
#     @assert mod(ntimesteps, nsave) == 0

#     if nwrite > 0
#         @assert mod(nwrite, nsave) == 0
#         @assert mod(ntimesteps, nwrite) == 0
#     end

#     N = nsamples(problem) > 1 ? 2 : 1
#     nd = length(vec(problem.ics.q))
#     nm = length(vec(problem.ics.λ))
#     ni = nsamples(problem)
#     nt = div(ntimesteps, nsave)
#     nt = (nwrite == 0 ? nt : div(nwrite, nsave))
#     nw = (nwrite == 0 ? ntimesteps : nwrite)

#     @assert nd > 0
#     @assert nm > 0
#     @assert ni > 0
#     @assert nt ≥ 0
#     @assert nw ≥ 0

#     t = TimeSeries{TT}(nt, timestep(problem), nsave)
#     q = DataSeries(problem.ics.q, nt, ni)
#     λ = DataSeries(problem.ics.λ, nt, ni)
#     s = SolutionDAE{AT,TT,N}(nd, nm, nt, ni, t, q, λ, ntimesteps, nsave, nw, periodicity(problem))
#     set_initial_conditions!(s, problem)

#     return s
# end

# function SolutionDAE(t::TimeSeries{TT}, q::DataSeries{DT,N}, λ::DataSeries{DT,N}, ntimesteps::Int) where {DT,TT,N}
#     @assert ndims(q) >= ndims(λ)
#     @assert ntime(q) == ntime(λ)
#     @assert nsamples(q) == nsamples(λ)

#     # extract parameters
#     nd = length(q[begin])
#     nm = length(λ[begin])
#     ni = nsamples(q)
#     nt = t.n
#     ns = div(ntimesteps, nt)

#     @assert mod(ntimesteps, nt) == 0

#     # create solution
#     SolutionDAE{DT,TT,N}(nd, nm, nt, ni, t, q, λ, ntimesteps, ns, 0)
# end

# function SolutionDAE(h5io::SolutionHDF5)
#     # get hdf5 file descriptor
#     h5 = hdf5(h5io)

#     # read attributes
#     ntimesteps = read(attributes(h5)["ntime"])
#     nsave = read(attributes(h5)["nsave"])
#     nsamples = read(attributes(h5)["nsamples"])

#     # reading data arrays
#     t = TimeSeries(read(h5["t"]), nsave)
#     q = fromarray(DataSeries, read(h5["q"]), nsamples)
#     λ = fromarray(DataSeries, read(h5["λ"]), nsamples)

#     # create solution
#     SolutionDAE(t, q, λ, ntimesteps)
# end

# function SolutionDAE(file::String)
#     load(file) do h5io
#         SolutionDAE(h5io::SolutionHDF5)
#     end
# end


# Base.:(==)(sol1::SolutionDAE{DT1,TT1,N1}, sol2::SolutionDAE{DT2,TT2,N2}) where {DT1,TT1,N1,DT2,TT2,N2} = (
#     DT1 == DT2
#     && TT1 == TT2
#     && N1 == N2
#     && sol1.nd == sol2.nd
#     && sol1.nm == sol2.nm
#     && sol1.nt == sol2.nt
#     && sol1.ni == sol2.ni
#     && sol1.t == sol2.t
#     && sol1.q == sol2.q
#     && sol1.λ == sol2.λ
#     && sol1.ntime == sol2.ntime
#     && sol1.nsave == sol2.nsave
#     && sol1.nwrite == sol2.nwrite
#     && sol1.offset == sol2.offset
#     && sol1.counter == sol2.counter
#     && sol1.periodicity == sol2.periodicity)

# @inline GeometricSolutions.counter(sol::SolutionDAE) = sol.counter
# @inline GeometricSolutions.offset(sol::SolutionDAE) = sol.offset
# @inline GeometricSolutions.lastentry(sol::SolutionDAE) = sol.ni == 1 ? sol.counter[1] - 1 : sol.counter .- 1

# @inline GeometricBase.nsamples(sol::SolutionDAE) = sol.ni
# @inline GeometricBase.nsave(sol::SolutionDAE) = sol.nsave
# @inline GeometricBase.ntime(sol::SolutionDAE) = sol.ntime
# @inline GeometricBase.timesteps(sol::SolutionDAE) = sol.t
# @inline GeometricBase.eachtimestep(sol::SolutionDAE) = 1:sol.nt*sol.nsave
# @inline GeometricBase.periodicity(sol::SolutionDAE) = sol.periodicity


# function get_initial_conditions!(sol::SolutionDAE{AT,TT}, asol::AtomicSolutionDAE{DT,TT,AT}, k, n = 1) where {DT,TT,AT<:AbstractArray{DT}} # TODO
# function get_initial_conditions!(sol::SolutionDAE, asol::AtomicSolutionDAE, k, n = 1)
#     get_solution!(sol, asol.q, asol.λ, n - 1, k)
#     asol.t = sol.t[n-1]
#     asol.q̃ .= 0
# end

# function get_initial_conditions!(sol::SolutionDAE{DT}, q::DT, λ::DT, k, n = 1) where {DT}
#     get_data!(sol.q, q, n - 1, k)
#     get_data!(sol.λ, λ, n - 1, k)
# end

# function get_initial_conditions(sol::SolutionDAE, k, n = 1)
#     get_solution(sol, n - 1, k)
# end

# function set_solution!(sol::SolutionDAE, t, q, λ, n, k = 1)
#     set_solution!(sol, q, λ, n, k)
# end

# function set_solution!(sol::SolutionDAE{AT,TT}, asol::AtomicSolutionDAE{DT,TT,AT}, n, k = 1) where {DT,TT,AT<:AbstractArray{DT}}
#     set_solution!(sol, asol.t, asol.q, asol.λ, n, k)
# end

# function set_solution!(sol::SolutionDAE{AT}, q::AT, λ::AT, n, k = 1) where {AT}
#     @assert n <= sol.ntime
#     @assert k <= sol.ni
#     if mod(n, sol.nsave) == 0
#         if sol.counter[k] > sol.nt
#             @error("Solution overflow. Call write_to_hdf5() and reset!() before continuing the simulation.")
#         end
#         set_data!(sol.q, q, sol.counter[k], k)
#         set_data!(sol.λ, λ, sol.counter[k], k)
#         sol.counter[k] += 1
#     end
# end

# function get_solution!(sol::SolutionDAE{AT}, q::AT, λ::AT, n, k = 1) where {AT} # TODO
# function get_solution!(sol::SolutionDAE, q, λ, n, k = 1)
#     q .= sol.q[n]
#     λ .= sol.λ[n]
# end

# function get_solution(sol::SolutionDAE, n, k = 1)
#     @assert k == 1
#     (sol.t[n], sol.q[n], sol.λ[n])
# end

# function get_solution(sol::SolutionDAE{AT,TT,1}, n, k = 1) where {AT,TT}
#     @assert k == 1
#     (sol.t[n], sol.q[n], sol.λ[n])
# end

# function get_solution(sol::SolutionDAE{AT,TT,2}, n, k = 1) where {AT,TT}
#     (sol.t[n], sol.q[n, k], sol.λ[n, k])
# end

# function GeometricBase.reset!(sol::SolutionDAE)
#     reset!(sol.q)
#     reset!(sol.λ)
#     compute_timeseries!(sol.t, sol.t[end])
#     sol.counter .= 1
#     sol.offset += sol.nt
# end

# function Base.close(solution::SolutionDAE)
#     close(solution.h5)
# end
