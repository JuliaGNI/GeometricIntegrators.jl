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
abstract type SolutionDAE{dType, tType, N} <: DeterministicSolution{dType, tType, N} end

# Create SolutionPDAEs with serial and parallel data structures.
for (TSolution, TDataSeries, Tdocstring) in
    ((:SSolutionDAE, :SDataSeries, "Serial Solution of a differential algebraic equation."),
     (:PSolutionDAE, :PDataSeries, "Parallel Solution of a differential algebraic equation."))
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
            periodicity::Vector{dType}
            h5::HDF5File

            function $TSolution{dType, tType, N}(nd, nm, nt, ni, t, q, λ, ntime, nsave, nwrite, periodicity=zeros(dType, nd)) where {dType <: Number, tType <: Real, N}
                new(nd, nm, nt, ni, t, q, λ, ntime, nsave, nwrite, zeros(Int, ni), 0, periodicity)
            end
        end

        function $TSolution(equation::DAE{DT,TT}, Δt::TT, ntime::Int; nsave::Int=DEFAULT_NSAVE, nwrite::Int=DEFAULT_NWRITE, filename=nothing) where {DT,TT}
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
            @assert nt ≥ 0
            @assert nw ≥ 0

            t = TimeSeries{TT}(nt, Δt, nsave)
            q = $TDataSeries(DT, nd, nt, ni)
            λ = $TDataSeries(DT, nm, nt, ni)
            s = $TSolution{DT,TT,N}(nd, nm, nt, ni, t, q, λ, ntime, nsave, nw, periodicity(equation))
            set_initial_conditions!(s, equation)

            if !isnothing(filename)
                isfile(filename) ? @warn("Overwriting existing HDF5 file.") : nothing
                create_hdf5!(s, filename)
            end

            return s
        end

        function $TSolution(t::TimeSeries{TT}, q::$TDataSeries{DT,N}, λ::$TDataSeries{DT,N}, ntime::Int) where {DT,TT,N}
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
            $TSolution{DT,TT,N}(nd, nm, nt, ni, t, q, λ, ntime, ns, 0)
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
            λ = $TDataSeries(read(h5["λ"]))

            # need to close the file
            close(h5)

            # create solution
            $TSolution(t, q, λ, ntime)
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
@inline timesteps(sol::SolutionDAE)  = sol.t
@inline ntime(sol::SolutionDAE) = sol.ntime
@inline nsave(sol::SolutionDAE) = sol.nsave
@inline offset(sol::SolutionDAE) = sol.woffset
@inline CommonFunctions.periodicity(sol::SolutionDAE) = sol.periodicity


"Create AtomicSolution for DAE."
function AtomicSolution(solution::SolutionDAE{DT,TT}) where {DT,TT}
    AtomicSolutionDAE(DT, TT, solution.nd, solution.nm)
end


function set_initial_conditions!(sol::SolutionDAE, equ::DAE)
    set_initial_conditions!(sol, equ.t₀, equ.q₀, equ.λ₀)
end

function set_initial_conditions!(sol::SolutionDAE{DT,TT}, t₀::TT, q₀::Array{DT}, λ₀::Array{DT}) where {DT,TT}
    set_data!(sol.q, q₀, 0)
    set_data!(sol.λ, λ₀, 0)
    compute_timeseries!(sol.t, t₀)
    sol.counter .= 1
end

function get_initial_conditions!(sol::SolutionDAE{DT,TT}, asol::AtomicSolutionDAE{DT,TT}, k, n=1) where {DT,TT}
    get_solution!(sol, asol.q, asol.λ, n-1, k)
    asol.t  = sol.t[n-1]
    asol.q̃ .= 0
end

function get_initial_conditions!(sol::SolutionDAE{DT,TT}, q::Vector{DT}, λ::Vector{DT}, k, n=1) where {DT,TT}
    get_data!(sol.q, q, n-1, k)
    get_data!(sol.λ, λ, n-1, k)
end

function get_initial_conditions(sol::SolutionDAE, k, n=1)
    get_solution(sol, n-1, k)
end

function set_solution!(sol::SolutionDAE, t, q, λ, n, k=1)
    set_solution!(sol, q, λ, n, k)
end

function set_solution!(sol::SolutionDAE{DT,TT}, asol::AtomicSolutionDAE{DT,TT}, n, k=1) where {DT,TT}
    set_solution!(sol, asol.t, asol.q, asol.λ, n, k)
end

function set_solution!(sol::SolutionDAE{DT,TT}, q::Vector{DT}, λ::Vector{DT}, n, k=1) where {DT,TT}
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

function get_solution!(sol::SolutionDAE{DT,TT}, q::SolutionVector{DT}, λ::SolutionVector{DT}, n, k=1) where {DT,TT}
    for i in eachindex(q) q[i] = sol.q[i, n, k] end
    for i in eachindex(λ) λ[i] = sol.λ[i, n, k] end
end

function get_solution(sol::SolutionDAE, n, k=1)
    (sol.t[n], sol.q[:, n, k], sol.λ[:, n, k])
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
