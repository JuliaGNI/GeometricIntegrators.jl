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

        @assert CONV==:strong || (CONV==:weak && K==0)
        @assert ntime==nt*nsave
        @assert q.ni == p.ni
        @assert q.nd == p.nd
        @assert q.nt == p.nt == t.n
        @assert W.ns == q.ni

        new{DT,TT,NQ,NW,CONV}(nd, nm, nt, ns, t, q, p, W, K, ntime, nsave, zeros(Int, ns), 0)
    end

    function SolutionPSDE(nd::Int, nm::Int, nt::Int, ns::Int, ni::Int, Δt::tType,
                W::WienerProcess{dType,tType,NW,CONV}, K::Int, ntime::Int, nsave::Int) where {dType <: Number, tType <: Real, NW, CONV}

        @assert CONV==:strong || (CONV==:weak && K==0)
        @assert nd > 0
        @assert ns > 0
        @assert ni > 0
        @assert ni == 1 || ns == 1
        @assert nsave > 0
        @assert ntime == 0 || ntime ≥ nsave
        @assert mod(ntime, nsave) == 0

        @assert NW ∈ (2,3)

        t = TimeSeries{tType}(nt, Δt, nsave)
        q = SDataSeries(dType, nd, nt, max(ns,ni))
        p = SDataSeries(dType, nd, nt, max(ns,ni))
        NQ = ns==ni==1 ? 2 : 3

        new{dType, tType, NQ, NW, CONV}(nd, nm, nt, max(ns,ni), t, q, p, W, K, ntime, nsave, zeros(Int, max(ns,ni)), 0)
    end

end


function SolutionPSDE(equation::Union{PSDE{DT,TT},SPSDE{DT,TT}}, Δt::TT, ntime::Int, nsave::Int=1; K::Int=0, conv=DEFAULT_SCONV) where {DT,TT}
    nd = equation.d
    nm = equation.m
    ns = equation.ns
    ni = equation.ni
    nt = div(ntime, nsave)

    # Holds the Wiener process data for ALL computed time steps
    # Wiener process increments are automatically generated here
    W = WienerProcess(DT, nm, ntime, ns, Δt, conv)

    s = SolutionPSDE(nd, nm, nt, ns, ni, Δt, W, K, ntime, nsave)
    set_initial_conditions!(s, equation)

    return s
end


function SolutionPSDE(equation::Union{PSDE{DT,TT},SPSDE{DT,TT}}, Δt::TT, dW::Array{DT, NW}, dZ::Array{DT, NW}, ntime::Int, nsave::Int=1; K::Int=0, conv=DEFAULT_SCONV) where {DT,TT,NW}
    nd = equation.d
    nm = equation.m
    ns = equation.ns
    ni = equation.ni
    nt = div(ntime, nsave)

    @assert size(dW) == size(dZ)
    @assert NW ∈ (2,3)

    if NW==2
        @assert nm==size(dW,1)
        @assert ntime==size(dW,2)
        @assert ns==1
    elseif NW==3
        @assert nm==size(dW,1)
        @assert ntime==size(dW,2)
        @assert ns==size(dW,3)
    end

    # Holds the Wiener process data for ALL computed time steps
    # Wiener process increments are prescribed by the arrays dW and dZ
    W = WienerProcess(Δt, dW, dZ, conv)

    s = SolutionPSDE(nd, nm, nt, ns, ni, Δt, W, K, ntime, nsave)
    set_initial_conditions!(s, equation)

    return s
end




# If the Wiener process W data are not available, creates a one-element zero array instead
# For instance used when reading a file with no Wiener process data saved
function SolutionPSDE(t::TimeSeries{TT}, q::SDataSeries{DT,NQ}, p::SDataSeries{DT,NQ}; K::Int=0, conv=DEFAULT_SCONV) where {DT,TT,NQ}
    # extract parameters
    nd = q.nd
    ns = q.ni
    nt = t.n
    nsave = t.step
    ntime = nt*nsave

    W = WienerProcess(t.Δt, [0.0], [0.0], conv)

    # create solution
    SolutionPSDE(nd, W.nd, nt, ns, 1, t, q, p, W, K, ntime, nsave)
end


function SolutionPSDE(file::String)
    # open HDF5 file
    get_config(:verbosity) > 1 ? @info("Reading HDF5 file ", file) : nothing
    h5 = h5open(file, "r")

    # read attributes
    nsave = read(attrs(h5)["nsave"])
    ni    = read(attrs(h5)["ni"])

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

    if ndims(q_array)==3 && ni>1
        q = SStochasticDataSeries(q_array,IC=true)
        p = SStochasticDataSeries(p_array,IC=true)
    else
        q = SStochasticDataSeries(q_array)
        p = SStochasticDataSeries(p_array)
    end

    # create solution
    if W_exists == true
        SolutionPSDE(t, q, p, W, K=K)
    else
        SolutionPSDE(t, q, p, K=K, conv=conv)
    end

end


timesteps(sol::SolutionPSDE)  = sol.t.t
ntime(sol::SolutionPSDE) = sol.ntime
nsave(sol::SolutionPSDE) = sol.nsave
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
    get_solution!(sol, asol.q, n-1, k)
    get_solution!(sol, asol.p, n-1, k)
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


function get_solution!(sol::SolutionPSDE{DT,TT}, q::SolutionVector{DT}, p::SolutionVector{DT}, n, k) where {DT,TT}
    for i in eachindex(q) q[i] = sol.q[i, n, k] end
    for i in eachindex(p) p[i] = sol.p[i, n, k] end
end

function get_solution(sol::SolutionPSDE, n, k)
    (sol.t[n], sol.q[:, n, k], sol.p[:, n, k])
end

function set_solution!(sol::SolutionPSDE, t, q, p, n, k)
    set_solution!(sol, q, p, n, k)
end

function set_solution!(sol::SolutionPSDE{DT,TT}, asol::AtomicSolutionPSDE{DT,TT}, n, k) where {DT,TT}
    set_solution!(sol, asol.t, asol.q, asol.p, n, k)
end

function set_solution!(sol::SolutionPSDE{DT,TT}, q::SolutionVector{DT}, p::SolutionVector{DT}, n, k) where {DT,TT}
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


function CommonFunctions.reset!(sol::SolutionPSDE)
    reset!(sol.q)
    reset!(sol.p)
    compute_timeseries!(sol.t, sol.t[end])
    generate_wienerprocess!(sol.W)
    sol.counter .= 1
    sol.woffset += sol.nt
end

function Base.close(solution::SolutionPSDE)
    # close(solution.t)
    # close(solution.q)
    close(solution.h5)
end



"""
Creates HDF5 file and initialises datasets for SDE solution object.
  It is implemented as one fucntion for all NQ and NW cases, rather than several
  separate cases as was done for SolutionODE.
  nt - the total number of time steps to store
  ntime - the total number of timesteps to be computed
"""
function create_hdf5(solution::SolutionPSDE{DT,TT,NQ,NW}, file::AbstractString,  nt::Int=solution.nt, ntime::Int=solution.ntime; save_W=true) where {DT,TT,NQ,NW}
    @assert nt ≥ 1
    @assert ntime ≥ 1

    # create HDF5 file
    h5 = createHDF5(solution, file)

    # Adding the attributes specific to SolutionSDE that were not added above
    save_attributes(solution, h5)
    attrs(h5)["nt"] = nt
    attrs(h5)["ntime"] = ntime

    # create dataset
    # nt and ntime can be used to set the expected total number of timesteps to be saved,
    # so that the size of the array does not need to be adapted dynamically.
    # Right now, it has to be set as dynamical size adaptation is not yet
    # working. The default value is the size of the solution structure.
    if NQ==2
        q = d_create(h5, "q", datatype(DT), dataspace(solution.nd, nt+1), "chunk", (solution.nd,1))
        p = d_create(h5, "p", datatype(DT), dataspace(solution.nd, nt+1), "chunk", (solution.nd,1))
        # copy initial conditions
        q[:, 1] = solution.q[:, 0]
        p[:, 1] = solution.p[:, 0]
    elseif NQ==3
        if solution.ns>1
            q = d_create(h5, "q", datatype(DT), dataspace(solution.nd, nt+1, solution.ns), "chunk", (solution.nd,1,1))
            p = d_create(h5, "p", datatype(DT), dataspace(solution.nd, nt+1, solution.ns), "chunk", (solution.nd,1,1))
        else
            q = d_create(h5, "q", datatype(DT), dataspace(solution.nd, nt+1, solution.ni), "chunk", (solution.nd,1,1))
            p = d_create(h5, "p", datatype(DT), dataspace(solution.nd, nt+1, solution.ni), "chunk", (solution.nd,1,1))
        end
        # copy initial conditions
        q[:, 1, :] = solution.q[:, 0, :]
        p[:, 1, :] = solution.p[:, 0, :]
    end

    if save_W==true
        # creating datasets to store the Wiener process increments
        if NW==2
            dW = d_create(h5, "ΔW", datatype(DT), dataspace(solution.nm, ntime), "chunk", (solution.nm,1))
            dZ = d_create(h5, "ΔZ", datatype(DT), dataspace(solution.nm, ntime), "chunk", (solution.nm,1))
        elseif NW==3
            dW = d_create(h5, "ΔW", datatype(DT), dataspace(solution.nm, ntime, solution.ns), "chunk", (solution.nm,1,1))
            dZ = d_create(h5, "ΔZ", datatype(DT), dataspace(solution.nm, ntime, solution.ns), "chunk", (solution.nm,1,1))
        end
    end

    # Creating a dataset for storing the time series
    t = d_create(h5, "t", datatype(TT), dataspace((nt+1,)), "chunk", (1,))
    t[1] = solution.t[0]

    return h5
end


"""
Append solution to HDF5 file.
  offset - start writing q at the position offset+2
  offset2- start writing ΔW, ΔZ at the position offset2+1
"""
function CommonFunctions.write_to_hdf5(solution::SolutionPSDE{DT,TT,NQ,NW}, h5::HDF5File, offset=0, offset2=offset) where {DT,TT,NQ,NW}
    # set convenience variables and compute ranges
    d   = solution.nd
    m   = solution.nm
    n   = solution.nt
    s   = solution.ns
    i   = solution.ni
    ntime = solution.ntime
    j1  = offset+2
    j2  = offset+1+n
    jw1 = offset2+1
    jw2 = offset2+ntime

    # # extend dataset if necessary
    # if size(x, 2) < j2
    #     set_dims!(x, (d, j2))
    # end

    # saving the time time series
    h5["t"][j1:j2] = solution.t[1:n]

    # copy data from solution to HDF5 dataset
    if NQ==2
        h5["q"][:, j1:j2] = solution.q[:, 1:n]
        h5["p"][:, j1:j2] = solution.p[:, 1:n]
    elseif NQ==3
        h5["q"][:, j1:j2, :] = solution.q[:, 1:n, :]
        h5["p"][:, j1:j2, :] = solution.p[:, 1:n, :]
    end

    if exists(h5, "ΔW") && exists(h5, "ΔZ")
        # copy the Wiener process increments from solution to HDF5 dataset
        if NW==2
            h5["ΔW"][:,jw1:jw2] = solution.W.ΔW.d[:,1:ntime]
            h5["ΔZ"][:,jw1:jw2] = solution.W.ΔZ.d[:,1:ntime]
        elseif NW==3
            h5["ΔW"][:,jw1:jw2,:] = solution.W.ΔW.d[:,1:ntime,:]
            h5["ΔZ"][:,jw1:jw2,:] = solution.W.ΔZ.d[:,1:ntime,:]
        end
    end


    return nothing
end
