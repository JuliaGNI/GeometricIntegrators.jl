"""
`SolutionSDE`: Solution of a stochastic differential equation

Contains all fields necessary to store the solution of an SDE.

### Fields

* `nd`: dimension of the dynamical variable ``q``
* `nm`: dimension of the Wiener process
* `nt`: number of time steps to store
* `ns`: number of sample paths
* `ni`: number of initial conditions
* `t`:  time steps
* `q`:  solution `q[nd, nt+1, ns, ni]` with `q[:,0,:,:]` the initial conditions
* `W`:  Wiener process driving the stochastic process q
* `ntime`: number of time steps to compute
* `nsave`: save every nsave'th time step

"""
mutable struct SolutionSDE{dType, tType, NQ, NW} <: StochasticSolution{dType, tType, NQ, NW}
    nd::Int
    nm::Int
    nt::Int
    ns::Int
    ni::Int
    t::TimeSeries{tType}
    q::SStochasticDataSeries{dType,NQ}
    W::WienerProcess{dType,tType,NW}
    ntime::Int
    nsave::Int
    counter::Int
end


function SolutionSDE(equation::SDE{DT,TT,VT,BT}, Δt::TT, ntime::Int, nsave::Int=1) where {DT,TT,VT,BT}
    nd = equation.d
    nm = equation.m
    ns = equation.ns
    ni = equation.n
    nt = div(ntime, nsave)

    if nd==ns==ni==1
        NQ = 1
    elseif ns==ni==1
        NQ = 2
    elseif ns==1 || ni==1
        NQ = 3
    else
        NQ = 4
    end

    @assert DT <: Number
    @assert TT <: Real
    @assert nd > 0
    @assert ns > 0
    @assert ni > 0
    @assert nsave > 0
    @assert ntime == 0 || ntime ≥ nsave
    @assert mod(ntime, nsave) == 0

    t = TimeSeries{TT}(nt, Δt, nsave)

    q = SStochasticDataSeries(DT, nd, nt, ns, ni)

    # Holds the Wiener process data for ALL computed time steps
    # Wiener process increments are automatically generated here
    W = WienerProcess(DT, nm, ntime, ns, Δt)
    NW = ndims(W.ΔW)

    s = SolutionSDE{DT,TT,NQ,NW}(nd, nm, nt, ns, ni, t, q, W, ntime, nsave, 0)
    set_initial_conditions!(s, equation)
    return s
end


function SolutionSDE(equation::SDE{DT,TT,VT,BT}, Δt::TT, dW::Array{DT, NW}, dZ::Array{DT, NW}, ntime::Int, nsave::Int=1) where {DT,TT,VT,BT,NW}
    nd = equation.d
    nm = equation.m
    ns = equation.ns
    ni = equation.n
    nt = div(ntime, nsave)

    if nd==ns==ni==1
        NQ = 1
    elseif ns==ni==1
        NQ = 2
    elseif ns==1 || ni==1
        NQ = 3
    else
        NQ = 4
    end

    @assert size(dW) == size(dZ)

    if NW==1
        @assert ns==nm==1
        @assert ntime==length(dW)
    elseif NW==2
        @assert nm==size(dW,1)
        @assert ntime==size(dW,2)
        @assert ns==1
    elseif NW==3
        @assert nm==size(dW,1)
        @assert ntime==size(dW,2)
        @assert ns==size(dW,3)
    end

    @assert DT <: Number
    @assert TT <: Real
    @assert nd > 0
    @assert ns > 0
    @assert ni > 0
    @assert nsave > 0
    @assert ntime == 0 || ntime ≥ nsave
    @assert mod(ntime, nsave) == 0


    t = TimeSeries{TT}(nt, Δt, nsave)

    q = SStochasticDataSeries(DT, nd, nt, ns, ni)

    # Holds the Wiener process data for ALL computed time steps
    # Wiener process increments are prescribed by the arrays ΔW and ΔZ
    W = WienerProcess(Δt, dW, dZ)

    s = SolutionSDE{DT,TT,NQ,NW}(nd, nm, nt, ns, ni, t, q, W, ntime, nsave, 0)
    set_initial_conditions!(s, equation)
    return s
end


function SolutionSDE(t::TimeSeries{TT}, q::SStochasticDataSeries{DT,NQ}, W::WienerProcess{DT,TT,NW}) where {DT,TT,NQ,NW}
    # extract parameters
    nd = q.nd
    ns = q.ns
    ni = q.ni
    nt = t.n
    nm = W.nd
    ntime = W.nt
    nsave = t.step

    @assert ntime==nt*nsave
    @assert q.nt == nt
    @assert W.ns == q.ns

    # create solution
    SolutionSDE{DT,TT,NQ,NW}(nd, nm, nt, ns, ni, t, q, W, ntime, nsave, 0)
end


function SolutionSDE(file::String)
    # open HDF5 file
    info("Reading HDF5 file ", file)
    h5 = h5open(file, "r")

    # read attributes
    nsave = read(attrs(h5)["nsave"])
    ni    = read(attrs(h5)["ni"])

    # reading data arrays
    t = TimeSeries(read(h5["t"]), nsave)
    W = WienerProcess(t.Δt, read(h5["ΔW"]), read(h5["ΔZ"]))

    q_array = read(h5["q"])

    close(h5)

    if ndims(q_array)==3 && ni>1
        q = SStochasticDataSeries(q_array,IC=true)
    else
        q = SStochasticDataSeries(q_array)
    end

    # create solution
    SolutionSDE(t, q, W)
end


time(sol::SolutionSDE)  = sol.t.t
ntime(sol::SolutionSDE) = sol.ntime
nsave(sol::SolutionSDE) = sol.nsave


function set_initial_conditions!(sol::SolutionSDE{DT,TT}, equ::SDE{DT,TT}) where {DT,TT}
    set_initial_conditions!(sol, equ.t₀, equ.q₀)
end


function set_initial_conditions!(sol::SolutionSDE{DT,TT}, t₀::TT, q₀::Union{Array{DT}, Array{Double{DT}}}) where {DT,TT}
    # Sets the initial conditions sol.q[0] with the data from q₀
    # q₀ may be 1D (nd elements - single deterministic initial condition),
    # 2D (nd x ns or nd x ni matrix - single random or multiple deterministic initial condition),
    # or 3D (nd x ns x ni matrix - multiple random initial condition)
    set_data!(sol.q, q₀, 0)
    compute_timeseries!(sol.t, t₀)
end


# copies the m-th initial condition for the k-th sample path from sol.q to q
function get_initial_conditions!(sol::SolutionSDE{DT,TT}, q::Union{Vector{DT}, Vector{Double{DT}}}, k, m) where {DT,TT}

    @assert k ≤ sol.ns
    @assert m ≤ sol.ni

    N = ndims(sol.q)

    if N==1
        # 1D space, 1 sample path and 1 initial condition, k==m==1
        q[1] = get_data!(sol.q, 0)
    elseif N==2
        # Multidimensional space, 1 sample path and 1 initial condition, k==m==1
        get_data!(sol.q, q, 0)
    elseif N==3
        # Multidimensional space, with either 1 sample path or 1 initial condition
        if sol.ns==1
            # 1 sample path, k==1, reading the m-th initial condition
            get_data!(sol.q, q, 0, m)
        else
            #1 initial condition, m==1, reading the k-th sample path
            get_data!(sol.q, q, 0, k)
        end
    elseif N==4
        # Multidimensional space, multiple sample paths and initial conditions
        get_data!(sol.q, q, 0, k, m)
    end

end


function copy_solution!(sol::SolutionSDE{DT,TT,NQ,NW}, q::Union{Vector{DT}, Vector{Double{DT}}}, n, k, m) where {DT,TT,NQ,NW}

    if mod(n, sol.nsave) == 0

        if NQ ∈ (1,2)
            # Single sample path and a single initial condition, k==m==1
            @assert k==m==1
            set_data!(sol.q, q, div(n, sol.nsave))
        elseif NQ==3
            if sol.ni==1
                #Single initial condition, multiple sample paths, m==1
                @assert m==1
                set_data!(sol.q, q, div(n, sol.nsave), k)
            else
                #Single sample path, multiple initial conditions, k==1
                @assert k==1
                set_data!(sol.q, q, div(n, sol.nsave), m)
            end
        elseif NQ==4
            # Multiple sample paths and initial conditions
            set_data!(sol.q, q, div(n, sol.nsave), k, m)
        end

        sol.counter += 1
    end
end


function reset!(sol::SolutionSDE)
    reset!(sol.q)
    compute_timeseries!(sol.t, sol.t[end])
    generate_wienerprocess!(sol.W)
    sol.counter = 0
end


# "Creates HDF5 file and initialises datasets for SDE solution object."
# It is implemented as one fucntion for all NQ and NW cases, rather than several
# separate cases as was done for SolutionODE.
function create_hdf5(solution::SolutionSDE{DT,TT,NQ,NW}, file::AbstractString, ntime::Int=1) where {DT,TT,NQ,NW}
    @assert ntime ≥ 1

    # create HDF5 file and save ntime, nsave as attributes, and t as the dataset called "t"
    h5 = createHDF5(solution, file)

    # Adding the attributes specific to SolutionSDE that were not added above
    attrs(h5)["nd"] = solution.nd
    attrs(h5)["nm"] = solution.nm
    attrs(h5)["ns"] = solution.ns
    attrs(h5)["ni"] = solution.ni
    attrs(h5)["nt"] = solution.nt

    # create dataset
    # ntime can be used to set the expected total number of timesteps
    # so that the size of the array does not need to be adapted dynamically.
    # Right now, it has to be set as dynamical size adaptation is not yet
    # working.
    if NQ==1
        # COULDN'T FIGURE OUT HOW TO USE d_create FOR A 1D ARRAY - the line below gives errors,
        # i.e., a 1x1 dataset is created, instead of an array
        #q = d_create(h5, "q", datatype(DT), dataspace(solution.nt+1,))
        #q[1] = solution.q.d[1]

        # INSTEAD, ALLOCATING A ZERO ARRAY q AND USING write TO CREATE A DATASET IN THE FILE
        q = zeros(DT,solution.nt+1)
        # copy initial conditions
        q[1] = solution.q.d[1]
        write(h5,"q",q)
    elseif NQ==2
        q = d_create(h5, "q", datatype(DT), dataspace(solution.nd, solution.nt+1), "chunk", (solution.nd,1))
        # copy initial conditions
        q[1:solution.nd, 1] = solution.q.d[1:solution.nd, 1]
    elseif NQ==3
        if solution.ns>1
            q = d_create(h5, "q", datatype(DT), dataspace(solution.nd, solution.nt+1, solution.ns), "chunk", (solution.nd,1,1))
        else
            q = d_create(h5, "q", datatype(DT), dataspace(solution.nd, solution.nt+1, solution.ni), "chunk", (solution.nd,1,1))
        end
        # copy initial conditions
        q[:, 1, :] = solution.q.d[:, 1, :]
    else
        q = d_create(h5, "q", datatype(DT), dataspace(solution.nd, solution.nt+1, solution.ns, solution.ni), "chunk", (solution.nd,1,1,1))
        # copy initial conditions
        q[:, 1, :, :] = solution.q.d[:, 1, :, :]
    end

    write(h5, "ΔW", solution.W.ΔW.d)
    write(h5, "ΔZ", solution.W.ΔZ.d)

    return h5
end


# "Append solution to HDF5 file."
function CommonFunctions.write_to_hdf5(solution::SolutionSDE{DT,TT,NQ,NW}, h5::HDF5.HDF5File, offset=0) where {DT,TT,NQ,NW}
    # set convenience variables and compute ranges
    d  = solution.nd
    n  = solution.nt
    s  = solution.ns
    i  = solution.ni
    j1 = offset+2
    j2 = offset+1+n

    # # extend dataset if necessary
    # if size(x, 2) < j2
    #     set_dims!(x, (d, j2))
    # end

    # copy data from solution to HDF5 dataset
    if NQ==1
        h5["q"][j1:j2] = solution.q.d[2:n+1]
    elseif NQ==2
        h5["q"][1:d, j1:j2] = solution.q.d[1:d, 2:n+1]
    elseif NQ==3
        h5["q"][:, j1:j2, :] = solution.q.d[:, 2:n+1,:]
    else
        h5["q"][:, j1:j2, :, :] = solution.q.d[:, 2:n+1, :, :]
    end

    return nothing
end
