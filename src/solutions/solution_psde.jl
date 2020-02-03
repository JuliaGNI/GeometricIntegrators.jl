"""
`SolutionPSDE`: Solution of a partitioned stochastic differential equation

Contains all fields necessary to store the solution of a PSDE or SPSDE

### Fields
* `conv`: type of the solution: :strong or :weak
* `nd`: dimension of the dynamical variable ``q``
* `nm`: dimension of the Wiener process
* `nt`: number of time steps to store
* `ns`: number of sample paths
* `ni`: number of initial conditions
* `t`:  time steps
* `q`:  solution `q[nd, nt+1, ns, ni]` with `q[:,0,:,:]` the initial conditions
* `p`:  solution `p[nd, nt+1, ns, ni]` with `p[:,0,:,:]` the initial conditions
* `W`:  Wiener process driving the stochastic processes q and p
* `K`:  integer parameter defining the truncation of the increments of the Wiener process (for strong solutions),
*       A = √(2 K Δt |log Δt|) due to Milstein & Tretyakov; if K=0 no truncation
* `ntime`: number of time steps to compute
* `nsave`: save every nsave'th time step

"""
mutable struct SolutionPSDE{dType, tType, NQ, NW, CONV} <: StochasticSolution{dType, tType, NQ, NW}
    nd::Int
    nm::Int
    nt::Int
    ns::Int
    t::TimeSeries{tType}
    q::SDataSeries{dType,NQ}
    p::SDataSeries{dType,NQ}
    W::WienerProcess{dType,tType,NW,CONV}
    K::Int
    ntime::Int
    nsave::Int
    nwrite::Int
    counter::Vector{Int}
    woffset::Int
    h5::HDF5File

    function SolutionPSDE(t::TimeSeries{TT}, q::SDataSeries{DT,NQ}, p::SDataSeries{DT,NQ}, W::WienerProcess{DT,TT,NW,CONV}; K::Int=0) where {DT,TT,NQ,NW,CONV}
        # extract parameters
        nd = q.nd
        ns = q.ni
        nt = q.nt
        nm = W.nd
        ntime = W.nt
        nsave = t.step
        nwrite = ntime

        @assert CONV==:strong || (CONV==:weak && K==0)
        @assert ntime==nt*nsave
        @assert q.ni == p.ni
        @assert q.nd == p.nd
        @assert q.nt == p.nt == t.n
        @assert W.ns == q.ni

        new{DT,TT,NQ,NW,CONV}(nd, nm, nt, ns, t, q, p, W, K, ntime, nsave, nwrite, zeros(Int, ns), 0)
    end

    function SolutionPSDE(nd::Int, nm::Int, nt::Int, ns::Int, ni::Int, Δt::tType,
                W::WienerProcess{dType,tType,NW,CONV}, K::Int, ntime::Int, nsave::Int, nwrite::Int) where {dType <: Number, tType <: Real, NW, CONV}

        @assert CONV==:strong || (CONV==:weak && K==0)
        @assert nd > 0
        @assert ns > 0
        @assert ni > 0
        @assert ni == 1 || ns == 1
        @assert nsave > 0
        @assert ntime == 0 || ntime ≥ nsave
        @assert mod(ntime, nsave) == 0

        if nwrite > 0
            @assert mod(nwrite, nsave) == 0
            @assert mod(ntime, nwrite) == 0
        end

        @assert NW ∈ (2,3)

        t = TimeSeries{tType}(nt, Δt, nsave)
        q = SDataSeries(dType, nd, nt, max(ns,ni))
        p = SDataSeries(dType, nd, nt, max(ns,ni))
        NQ = ns==ni==1 ? 2 : 3

        new{dType, tType, NQ, NW, CONV}(nd, nm, nt, max(ns,ni), t, q, p, W, K, ntime, nsave, nwrite, zeros(Int, max(ns,ni)), 0)
    end
end


function SolutionPSDE(equation::Union{PSDE{DT,TT},SPSDE{DT,TT}}, Δt::TT, ntime::Int, nsave::Int=DEFAULT_NSAVE, nwrite::Int=DEFAULT_NWRITE; K::Int=0, conv=DEFAULT_SCONV, filename=nothing) where {DT,TT}
    nd = equation.d
    nm = equation.m
    ns = equation.ns
    ni = equation.ni
    nt = div(ntime, nsave)
    nt = (nwrite == 0 ? nt : div(nwrite, nsave))
    nw = (nwrite == 0 ? ntime : nwrite)

    # Holds the Wiener process data for ALL computed time steps
    # Wiener process increments are automatically generated here
    W = WienerProcess(DT, nm, ntime, max(ni,ns), Δt, conv)

    s = SolutionPSDE(nd, nm, nt, ns, ni, Δt, W, K, ntime, nsave, nw)
    set_initial_conditions!(s, equation)

    if !isnothing(filename)
        isfile(filename) ? @warn("Overwriting existing HDF5 file.") : nothing
        create_hdf5(s, filename)
    end

    return s
end


function SolutionPSDE(equation::Union{PSDE{DT,TT},SPSDE{DT,TT}}, Δt::TT, dW::Array{DT, NW}, dZ::Array{DT, NW}, ntime::Int, nsave::Int=DEFAULT_NSAVE, nwrite::Int=DEFAULT_NWRITE; K::Int=0, conv=DEFAULT_SCONV, filename=nothing) where {DT,TT,NW}
    nd = equation.d
    nm = equation.m
    ns = equation.ns
    ni = equation.ni
    nt = div(ntime, nsave)
    nt = (nwrite == 0 ? nt : div(nwrite, nsave))
    nw = (nwrite == 0 ? ntime : nwrite)

    @assert size(dW) == size(dZ)
    @assert NW ∈ (2,3)

    @assert nm == size(dW,1)
    @assert ntime == size(dW,2)
    @assert ns == size(dW,3)

    # Holds the Wiener process data for ALL computed time steps
    # Wiener process increments are prescribed by the arrays dW and dZ
    W = WienerProcess(Δt, dW, dZ, conv)

    s = SolutionPSDE(nd, nm, nt, ns, ni, Δt, W, K, ntime, nsave, nw)
    set_initial_conditions!(s, equation)

    if !isnothing(filename)
        isfile(filename) ? @warn("Overwriting existing HDF5 file.") : nothing
        create_hdf5(s, filename)
    end

    return s
end




function SolutionPSDE(file::String)
    # open HDF5 file
    get_config(:verbosity) > 1 ? @info("Reading HDF5 file ", file) : nothing
    h5 = h5open(file, "r")

    # read attributes
    nsave = read(attrs(h5)["nsave"])
    ns    = read(attrs(h5)["ns"])

    # reading data arrays
    t = TimeSeries(read(h5["t"]), nsave)

    if exists(attrs(h5),"conv")
        conv = Symbol(read(attrs(h5)["conv"]))
    else
        conv = DEFAULT_SCONV
    end

    W_exists = exists(h5, "ΔW") && exists(h5, "ΔZ")

    if W_exists == true
        W = WienerProcess(t.Δt, read(h5["ΔW"]), read(h5["ΔZ"]), conv)
    end

    if exists(attrs(h5),"K")
        K = read(attrs(h5)["K"])
    else
        K=0
    end

    q_array = read(h5["q"])
    p_array = read(h5["p"])

    close(h5)

    q = SDataSeries(q_array)
    p = SDataSeries(p_array)

    # create solution
    if W_exists == true
        SolutionPSDE(t, q, p, W, K=K)
    else
        SolutionPSDE(t, q, p, K=K, conv=conv)
    end

end

Base.:(==)(sol1::SolutionPSDE{DT1,TT1,NQ1,NW1,C1}, sol2::SolutionPSDE{DT2,TT2,NQ2,NW2,C2}) where {DT1,TT1,NQ1,NW1,C1,DT2,TT2,NQ2,NW2,C2} = (
                                DT1 == DT2
                             && TT1 == TT2
                             && NQ1 == NQ2
                             && NW1 == NW2
                             && C1  == C2
                             && sol1.nd == sol2.nd
                             && sol1.nm == sol2.nm
                             && sol1.nt == sol2.nt
                             && sol1.ns == sol2.ns
                             && sol1.t  == sol2.t
                             && sol1.q  == sol2.q
                             && sol1.p  == sol2.p
                             && sol1.W  == sol2.W
                             && sol1.K  == sol2.K
                             && sol1.ntime == sol2.ntime
                             && sol1.nsave == sol2.nsave
                             && sol1.nwrite == sol2.nwrite
                             && sol1.counter == sol2.counter
                             && sol1.woffset == sol2.woffset)

hdf5(sol::SolutionPSDE) = sol.h5
timesteps(sol::SolutionPSDE)  = sol.t
ntime(sol::SolutionPSDE) = sol.ntime
nsave(sol::SolutionPSDE) = sol.nsave
offset(sol::SolutionPSDE) = sol.woffset
conv(sol::SolutionPSDE{DT,TT,NQ,NW,CONV}) where {DT,TT,NQ,NW,CONV} = CONV


function set_initial_conditions!(sol::SolutionPSDE, equ::Union{PSDE,SPSDE})
    set_initial_conditions!(sol, equ.t₀, equ.q₀, equ.p₀)
end


function set_initial_conditions!(sol::SolutionPSDE{DT,TT}, t₀::TT, q₀::AbstractArray{DT,1}, p₀::AbstractArray{DT,1}) where {DT,TT}
    # Sets the initial conditions sol.q[0] with the data from q₀
    # Here, q₀ is 1D (nd elements) representing a single deterministic or
    # multiple random initial conditions.
    # Similar for sol.p[0].
    if sol.ns == 1
        set_data!(sol.q, q₀, 0)
        set_data!(sol.p, p₀, 0)
    else
        for k in 1:sol.ns
            set_data!(sol.q, q₀, 0, k)
            set_data!(sol.p, p₀, 0, k)
        end
    end
    compute_timeseries!(sol.t, t₀)
    sol.counter .= 1
end

function set_initial_conditions!(sol::SolutionPSDE{DT,TT}, t₀::TT, q₀::AbstractArray{DT,2}, p₀::AbstractArray{DT,2}) where {DT,TT}
    # Sets the initial conditions sol.q[0] with the data from q₀
    # Here, q₀ is a 2D (nd x ni) matrix representing multiple deterministic initial conditions.
    # Similar for sol.p[0].
    @assert sol.ns == size(q₀,2) == size(p₀,2)
    set_data!(sol.q, q₀, 0)
    set_data!(sol.p, p₀, 0)
    compute_timeseries!(sol.t, t₀)
    sol.counter .= 1
end


function get_initial_conditions!(sol::SolutionPSDE{DT,TT}, asol::AtomicSolutionPSDE{DT,TT}, k, n=1) where {DT,TT}
    get_solution!(sol, asol.q, asol.p, n-1, k)
    asol.t  = sol.t[n-1]
    asol.q̃ .= 0
    asol.p̃ .= 0
end

function get_initial_conditions!(sol::SolutionPSDE{DT}, q::SolutionVector{DT}, p::SolutionVector{DT}, k, n=1) where {DT}
    get_solution!(sol, q, p, n-1, k)
end

function get_initial_conditions(sol::SolutionPSDE, k, n=1)
    get_solution(sol, n-1, k)
end


function get_solution!(sol::SolutionPSDE{DT,TT}, q::SolutionVector{DT}, p::SolutionVector{DT}, n, k=1) where {DT,TT}
    for i in eachindex(q) q[i] = sol.q[i, n, k] end
    for i in eachindex(p) p[i] = sol.p[i, n, k] end
end

function get_solution(sol::SolutionPSDE, n, k=1)
    (sol.t[n], sol.q[:, n, k], sol.p[:, n, k])
end

function set_solution!(sol::SolutionPSDE, t, q, p, n, k=1)
    set_solution!(sol, q, p, n, k)
end

function set_solution!(sol::SolutionPSDE{DT,TT}, asol::AtomicSolutionPSDE{DT,TT}, n, k=1) where {DT,TT}
    set_solution!(sol, asol.t, asol.q, asol.p, n, k)
end

function set_solution!(sol::SolutionPSDE{DT,TT}, q::SolutionVector{DT}, p::SolutionVector{DT}, n, k=1) where {DT,TT}
    @assert n <= sol.ntime
    @assert k <= sol.ns
    if mod(n, sol.nsave) == 0
        if sol.counter[k] > sol.nt
            @error("Solution overflow. Call write_to_hdf5() and reset!() before continuing the simulation.")
        end
        set_data!(sol.q, q, sol.counter[k], k)
        set_data!(sol.p, p, sol.counter[k], k)
        sol.counter[k] += 1
    end
end


# copy increments of the Brownian Process for multidimensional Brownian motion, 1 sample path
function get_increment(sol::SolutionPSDE{DT,TT,NQ,2}, k=1) where {DT,TT,NQ}
    @assert k==1
    return (sol.W.ΔW[:,n], sol.W.ΔZ[:,n-1])
end

# copy increments of the Brownian Process for multidimensional Brownian motion, r-th sample path
function get_increment(sol::SolutionPSDE{DT,TT,NQ,3}, k) where {DT,TT,NQ}
    return (sol.W.ΔW[:,n,k], sol.W.ΔZ[:,n-1,k])
end

# copy increments of the Brownian Process for multidimensional Brownian motion, 1 sample path
function get_increments!(sol::SolutionPSDE{DT,TT,NQ,2}, asol::AtomicSolutionPSDE{DT,TT}, n, k=1) where {DT,TT,NQ}
    @assert k==1
    for l = 1:sol.nm
        asol.ΔW[l] = sol.W.ΔW[l,n]
        asol.ΔZ[l] = sol.W.ΔZ[l,n]
    end
end

# copy increments of the Brownian Process for multidimensional Brownian motion, r-th sample path
function get_increments!(sol::SolutionPSDE{DT,TT,NQ,3}, asol::AtomicSolutionPSDE{DT,TT}, n, k) where {DT,TT,NQ}
    for l = 1:sol.nm
        asol.ΔW[l] = sol.W.ΔW[l,n,k]
        asol.ΔZ[l] = sol.W.ΔZ[l,n,k]
    end
end


function CommonFunctions.reset!(sol::SolutionPSDE)
    reset!(sol.q)
    reset!(sol.p)
    compute_timeseries!(sol.t, sol.t[end])
    generate_wienerprocess!(sol.W)
    sol.counter .= 1
    sol.woffset += sol.nt
end

function Base.close(solution::SolutionPSDE)
    close(solution.h5)
end
