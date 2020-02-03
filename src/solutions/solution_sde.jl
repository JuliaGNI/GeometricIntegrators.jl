"""
`SolutionSDE`: Solution of a stochastic differential equation

Contains all fields necessary to store the solution of an SDE.

### Fields
* `conv`: type of the solution: :strong or :weak
* `nd`: dimension of the dynamical variable ``q``
* `nm`: dimension of the Wiener process
* `nt`: number of time steps to store
* `ns`: number of sample paths
* `t`:  time steps
* `q`:  solution `q[nd, nt+1, ns]` with `q[:,0,:]` the initial conditions
* `W`:  Wiener process driving the stochastic process q
* `K`:  integer parameter defining the truncation of the increments of the Wiener process (for strong solutions),
*       A = √(2 K Δt |log Δt|) due to Milstein & Tretyakov; if K=0 no truncation
* `ntime`: number of time steps to compute
* `nsave`: save every nsave'th time step

"""
mutable struct SolutionSDE{dType, tType, NQ, NW, CONV} <: StochasticSolution{dType, tType, NQ, NW}
    nd::Int
    nm::Int
    nt::Int
    ns::Int
    t::TimeSeries{tType}
    q::SDataSeries{dType,NQ}
    W::WienerProcess{dType,tType,NW,CONV}
    K::Int
    ntime::Int
    nsave::Int
    nwrite::Int
    counter::Vector{Int}
    woffset::Int
    h5::HDF5File

    function SolutionSDE(t::TimeSeries{TT}, q::SDataSeries{DT,NQ}, W::WienerProcess{DT,TT,NW,CONV}; K::Int=0) where {DT,TT,NQ,NW,CONV}
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
        @assert q.nt == t.n
        @assert W.ns == q.ni

        new{DT,TT,NQ,NW,CONV}(nd, nm, nt, ns, t, q, W, K, ntime, nsave, nwrite, zeros(Int, ns), 0)
    end

    function SolutionSDE(nd::Int, nm::Int, nt::Int, ns::Int, ni::Int, Δt::tType,
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
        NQ = ns==ni==1 ? 2 : 3

        new{dType, tType, NQ, NW, CONV}(nd, nm, nt, max(ns,ni), t, q, W, K, ntime, nsave, nwrite, zeros(Int, max(ns,ni)), 0)
    end
end


function SolutionSDE(equation::SDE{DT,TT}, Δt::TT, ntime::Int, nsave::Int=DEFAULT_NSAVE, nwrite::Int=DEFAULT_NWRITE; K::Int=0, conv=DEFAULT_SCONV, filename=nothing) where {DT,TT}
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

    s = SolutionSDE(nd, nm, nt, ns, ni, Δt, W, K, ntime, nsave, nw)
    set_initial_conditions!(s, equation)

    if !isnothing(filename)
        isfile(filename) ? @warn("Overwriting existing HDF5 file.") : nothing
        create_hdf5(s, filename)
    end

    return s
end


function SolutionSDE(equation::SDE{DT,TT}, Δt::TT, dW::Array{DT, NW}, dZ::Array{DT, NW}, ntime::Int, nsave::Int=DEFAULT_NSAVE, nwrite::Int=DEFAULT_NWRITE; K::Int=0, conv=DEFAULT_SCONV, filename=nothing) where {DT,TT,NW}
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
    # Wiener process increments are prescribed by the arrays ΔW and ΔZ
    W = WienerProcess(Δt, dW, dZ, conv)

    s = SolutionSDE(nd, nm, nt, ns, ni, Δt, W, K, ntime, nsave, nw)
    set_initial_conditions!(s, equation)

    if !isnothing(filename)
        isfile(filename) ? @warn("Overwriting existing HDF5 file.") : nothing
        create_hdf5(s, filename)
    end

    return s
end


# function SolutionSDE(t::TimeSeries{TT}, q::SDataSeries{DT,NQ}, W::WienerProcess{DT,TT,NW,CONV}; K::Int=0) where {DT,TT,NQ,NW,CONV}
#     # extract parameters
#     nd = q.nd
#     ns = q.ni
#     nt = t.nt
#     nm = W.nd
#     ntime = W.nt
#     nsave = t.step
#
#     @assert conv==:strong || (conv==:weak && K==0)
#     @assert ntime==nt*nsave
#     @assert q.nt == nt
#     @assert W.ns == q.ns
#
#     # create solution
#     SolutionSDE{DT,TT,NQ,NW}(CONV, nd, nm, nt, ns, 1, t, q, W, K, ntime, nsave)
# end


# # If the Wiener process W data are not available, creates a one-element zero array instead
# # For instance used when reading a file with no Wiener process data saved
# function SolutionSDE(t::TimeSeries{TT}, q::SDataSeries{DT,NQ}; K::Int=0, conv=DEFAULT_SCONV) where {DT,TT,NQ}
#     # extract parameters
#     nd = q.nd
#     ni = q.ni
#     nt = q.nt
#     nsave = t.step
#     ntime = nt*nsave
#
#     W = WienerProcess(t.Δt, [0.0], [0.0], conv)
#
#     # create solution
#     SolutionSDE(nd, W.nd, nt, 1, ni, t, q, W, K, ntime, nsave, ntime)
# end


function SolutionSDE(file::String)
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

    close(h5)

    q = SDataSeries(q_array)

    # create solution
    if W_exists == true
        return SolutionSDE(t, q, W, K=K)
    else
        return SolutionSDE(t, q, K=K, conv=conv)
    end
end

Base.:(==)(sol1::SolutionSDE{DT1,TT1,NQ1,NW1,C1}, sol2::SolutionSDE{DT2,TT2,NQ2,NW2,C2}) where {DT1,TT1,NQ1,NW1,C1,DT2,TT2,NQ2,NW2,C2} = (
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
                             && sol1.W  == sol2.W
                             && sol1.K  == sol2.K
                             && sol1.ntime == sol2.ntime
                             && sol1.nsave == sol2.nsave
                             && sol1.nwrite == sol2.nwrite
                             && sol1.counter == sol2.counter
                             && sol1.woffset == sol2.woffset)

hdf5(sol::SolutionSDE) = sol.h5
timesteps(sol::SolutionSDE) = sol.t
ntime(sol::SolutionSDE) = sol.ntime
nsave(sol::SolutionSDE) = sol.nsave
offset(sol::SolutionSDE) = sol.woffset
conv(sol::SolutionSDE{DT,TT,NQ,NW,CONV}) where {DT,TT,NQ,NW,CONV} = CONV


function set_initial_conditions!(sol::SolutionSDE, equ::SDE)
    set_initial_conditions!(sol, equ.t₀, equ.q₀)
end


function set_initial_conditions!(sol::SolutionSDE{DT,TT}, t₀::TT, q₀::AbstractArray{DT,1}) where {DT,TT}
    # Sets the initial conditions sol.q[0] with the data from q₀
    # Here, q₀ is 1D (nd elements) representing a single deterministic or
    # multiple random initial conditions.
    if sol.ns == 1
        set_data!(sol.q, q₀, 0)
    else
        for k in 1:sol.ns
            set_data!(sol.q, q₀, 0, k)
        end
    end
    compute_timeseries!(sol.t, t₀)
    sol.counter .= 1
end


function set_initial_conditions!(sol::SolutionSDE{DT,TT}, t₀::TT, q₀::AbstractArray{DT,2}) where {DT,TT}
    # Sets the initial conditions sol.q[0] with the data from q₀
    # Here, q₀ is a 2D (nd x ni) matrix representing multiple deterministic initial conditions.
    @assert sol.ns == size(q₀,2)
    set_data!(sol.q, q₀, 0)
    compute_timeseries!(sol.t, t₀)
    sol.counter .= 1
end

function get_initial_conditions!(sol::SolutionSDE{DT,TT}, asol::AtomicSolutionSDE{DT,TT}, k, n=1) where {DT,TT}
    get_solution!(sol, asol.q, n-1, k)
    asol.t  = sol.t[n-1]
    asol.q̃ .= 0
end

# copies the m-th initial condition from sol.q to q
function get_initial_conditions!(sol::SolutionSDE{DT}, q::SolutionVector{DT}, k, n=1) where {DT}
    get_solution!(sol, q, n-1, k)
end

function get_initial_conditions(sol::SolutionSDE, k, n=1)
    get_solution(sol, n-1, k)
end

function get_solution!(sol::SolutionSDE{DT}, q::SolutionVector{DT}, n, k=1) where {DT}
    for i in eachindex(q) q[i] = sol.q[i, n, k] end
end

function get_solution(sol::SolutionSDE, n, k=1)
    (sol.t[n], sol.q[:, n, k])
end

function set_solution!(sol::SolutionSDE, t, q, n, k=1)
    set_solution!(sol, q, n, k)
end

function set_solution!(sol::SolutionSDE{DT,TT}, asol::AtomicSolutionSDE{DT,TT}, n, k=1) where {DT,TT}
    set_solution!(sol, asol.t, asol.q, n, k)
end

function set_solution!(sol::SolutionSDE{DT}, q::SolutionVector{DT}, n, k=1) where {DT}
    @assert n <= sol.ntime
    @assert k <= sol.ns
    if mod(n, sol.nsave) == 0
        if sol.counter[k] > sol.nt
            @error("Solution overflow. Call write_to_hdf5() and reset!() before continuing the simulation.")
        end
        set_data!(sol.q, q, sol.counter[k], k)
        sol.counter[k] += 1
    end
end

function CommonFunctions.reset!(sol::SolutionSDE)
    reset!(sol.q)
    compute_timeseries!(sol.t, sol.t[end])
    generate_wienerprocess!(sol.W)
    sol.counter .= 1
    sol.woffset += sol.nt
end

function Base.close(solution::SolutionSDE)
    # close(solution.t)
    # close(solution.q)
    close(solution.h5)
end


"""
Creates HDF5 file and initialises datasets for SDE solution object.
  It is implemented as one function for all NQ and NW cases, rather than several
  separate cases as was done for SolutionSDE.
  nt - the total number of time steps to store
  ntime - the total number of timesteps to be computed
"""
function create_hdf5(solution::SolutionSDE{DT,TT,NQ,NW}, file::AbstractString; save_W=true) where {DT,TT,NQ,NW}
    # create HDF5 file
    solution.h5 = createHDF5(solution, file)

    # Adding the attributes specific to SolutionSDE that were not added above
    save_attributes(solution, solution.h5)

    # create dataset
    # nt and ntime can be used to set the expected total number of timesteps to be saved,
    # so that the size of the array does not need to be adapted dynamically.
    # Right now, it has to be set as dynamical size adaptation is not yet
    # working. The default value is the size of the solution structure.
    if NQ==2
        q = d_create(solution.h5, "q", datatype(DT), dataspace(solution.nd, solution.nt+1), "chunk", (solution.nd,1))
        # copy initial conditions
        q[:,1] = solution.q[:,0]
    elseif NQ==3
        q = d_create(solution.h5, "q", datatype(DT), dataspace(solution.nd, solution.nt+1, solution.ns), "chunk", (solution.nd,1,1))
        # copy initial conditions
        q[:,1,:] = solution.q[:,0,:]
    end

    if save_W==true
        # creating datasets to store the Wiener process increments
        if NW==2
            dW = d_create(solution.h5, "ΔW", datatype(DT), dataspace(solution.nm, solution.ntime), "chunk", (solution.nm,1))
            dZ = d_create(solution.h5, "ΔZ", datatype(DT), dataspace(solution.nm, solution.ntime), "chunk", (solution.nm,1))
        elseif NW==3
            dW = d_create(solution.h5, "ΔW", datatype(DT), dataspace(solution.nm, solution.ntime, solution.ns), "chunk", (solution.nm,1,1))
            dZ = d_create(solution.h5, "ΔZ", datatype(DT), dataspace(solution.nm, solution.ntime, solution.ns), "chunk", (solution.nm,1,1))
        end
    end

    # Creating a dataset for storing the time series
    t = d_create(solution.h5, "t", datatype(TT), dataspace((solution.nt+1,)), "chunk", (1,))
    t[1] = solution.t[0]

    return solution.h5
end


function copy_solution_to_hdf5(solution::SolutionSDE{DT,TT,2,NW}, h5::HDF5File, j1, j2, n1, n2) where {DT,TT,NW}
    h5["q"][:, j1:j2] = solution.q[:, n1:n2]
end

function copy_solution_to_hdf5(solution::SolutionSDE{DT,TT,3,NW}, h5::HDF5File, j1, j2, n1, n2) where {DT,TT,NW}
    h5["q"][:, j1:j2, :] = solution.q[:, n1:n2, :]
end
