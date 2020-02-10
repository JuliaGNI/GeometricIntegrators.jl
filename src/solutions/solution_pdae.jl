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

        function $TSolution(equation::Union{IODE{DT,TT},VODE{DT,TT},PDAE{DT,TT},IDAE{DT,TT},VDAE{DT,TT}}, Δt::TT, ntime::Int; nsave::Int=DEFAULT_NSAVE, nwrite::Int=DEFAULT_NWRITE, filename=nothing) where {DT,TT}
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
            p = $TDataSeries(DT, nd, nt, ni)
            λ = $TDataSeries(DT, nm, nt, ni)
            s = $TSolution{DT,TT,N}(nd, nm, nt, ni, t, q, p, λ, ntime, nsave, nw)
            set_initial_conditions!(s, equation)

            if !isnothing(filename)
                isfile(filename) ? @warn("Overwriting existing HDF5 file.") : nothing
                create_hdf5!(s, filename)
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
                             && sol1.woffset == sol2.woffset)

hdf5(sol::SolutionPDAE)  = sol.h5
timesteps(sol::SolutionPDAE)  = sol.t
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

function get_solution!(sol::SolutionPDAE{DT,TT}, q::SolutionVector{DT}, p::SolutionVector{DT}, λ::SolutionVector{DT}, n, k=1) where {DT,TT}
    q .= sol.q[:, n, k]
    p .= sol.p[:, n, k]
    λ .= sol.λ[:, n, k]
end

function get_solution!(sol::SolutionPDAE{DT,TT}, q::SolutionVector{DT}, p::SolutionVector{DT}, n, k=1) where {DT,TT}
    q .= sol.q[:, n, k]
    p .= sol.p[:, n, k]
end

function get_solution(sol::SolutionPDAE, n, k=1)
    (sol.t[n], sol.q[:, n, k], sol.p[:, n, k], sol.λ[:, n, k])
end

function set_solution!(sol::SolutionPDAE, t::Real, q::Vector, p::Vector, n::Int, k::Int=1)
    set_solution!(sol, q, p, n, k)
end

function set_solution!(sol::SolutionPDAE, t::Real, q::Vector, p::Vector, λ::Vector, n::Int, k::Int=1)
    set_solution!(sol, q, p, λ, n, k)
end

function set_solution!(sol::SolutionPDAE{DT,TT}, asol::AtomicSolutionPDAE{DT,TT}, n::Int, k::Int=1) where {DT,TT}
    set_solution!(sol, asol.t, asol.q, asol.p, asol.λ, n, k)
end

function set_solution!(sol::SolutionPDAE{DT,TT}, q::SolutionVector{DT}, p::SolutionVector{DT}, λ::SolutionVector{DT}, n::Int, k::Int=1) where {DT,TT}
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

function set_solution!(sol::SolutionPDAE{DT,TT}, q::SolutionVector{DT}, p::SolutionVector{DT}, n::Int, k::Int=1) where {DT,TT}
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
    close(solution.h5)
end
