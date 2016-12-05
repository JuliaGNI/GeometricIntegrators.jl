"""
`SolutionODE`: Solution of an ordinary differential equation

Contains all fields necessary to store the solution of an ODE.

### Fields

* `nd`: dimension of the dynamical variable ``q``
* `nt`: number of time steps to store
* `n0`: number of initial conditions
* `t`:  time steps
* `x`:  solution `x[nd, nt+1, n0]` with `x[:,0,:]` the initial conditions
* `ntime`: number of time steps to compute
* `nsave`: save every nsave'th time step

"""
immutable SolutionODE{dType, tType} <: Solution{dType, tType, 3}
    nd::Int
    nt::Int
    n0::Int
    t::Timeseries{tType}
    x::Array{dType,3}
    ntime::Int
    nsave::Int

    function SolutionODE(nd, n0, ntime, nsave, Δt)
        @assert dType <: Number
        @assert tType <: Real
        @assert nd > 0
        @assert n0 > 0
        @assert nsave > 0
        @assert ntime ≥ nsave
        @assert mod(ntime, nsave) == 0

        nt = div(ntime, nsave)
        t = Timeseries{tType}(nt, Δt, nsave)
        x = zeros(dType, nd, nt+1, n0)
        new(nd, nt, n0, t, x, ntime, nsave)
    end
end

function SolutionODE{DT,TT,FT}(equation::ODE{DT,TT,FT}, Δt::TT, ntime::Int, nsave::Int=1)
    s = SolutionODE{DT,TT}(equation.d, equation.n, ntime, nsave, Δt)
    set_initial_conditions!(s, equation)
    return s
end

function set_initial_conditions!(solution::SolutionODE, equation::ODE)
    set_initial_conditions!(solution, equation.t₀, equation.q₀)
end

function set_initial_conditions!(solution::SolutionODE, t₀, q₀)
    simd_copy_yx_first_last!(q₀, solution, 0)
    solution.t[0] = t₀
    compute_timeseries!(solution.t)
end

function reset(s::SolutionODE)
    for k in 1:solution.n0
        for i in 1:solution.nd
            solution[i, 0, k] = solution[i, end, k]
        end
    end
    solution.t[0] = solution.t[end]
    compute_timeseries!(solution.t)
end

Base.indices(s::SolutionODE) = (1:s.nd, 0:s.nt, 1:s.n0)
Base.strides(s::SolutionODE) = (1, s.nd, s.nd*s.nt)
# Base.linearindexing{T<:SolutionODE}(::Type{T}) = LinearFast()

@inline function Base.getindex(s::SolutionODE, i::Int, j::Int, k::Int)
    @boundscheck checkbounds(s.x, i, j+1, k)
    @inbounds r = getindex(s.x, i, j+1, k)
    return r
end

@inline function Base.getindex(s::SolutionODE, i::Int, j::Int)
    if s.n0 == 1
        @boundscheck checkbounds(s.x, i, j+1, 1)
        @inbounds r = getindex(s.x, i, j+1, 1)
    else
        @boundscheck checkbounds(s.x, 1:s.nd, i+1, j)
        @inbounds r = getindex(s.x, 1:s.nd, i+1, j)
    end
    return r
end

@inline function Base.getindex(s::SolutionODE, k::Int)
    if s.n0 == 1
        @boundscheck checkbounds(s.x, 1:s.nd, k+1, 1)
        @inbounds r = getindex(s.x, 1:s.nd, k+1, 1)
    else
        @boundscheck checkbounds(s.x, 1:s.nd, 1:s.nt, k)
        @inbounds r = getindex(s.x, 1:s.nd, 1:s.nt, k)
    end
    return r
end

@inline function Base.setindex!(s::SolutionODE, x, i::Int, j::Int, k::Int)
    @boundscheck checkbounds(s.x, i, j+1, k)
    @inbounds setindex!(s.x, x, i, j+1, k)
end

@inline function Base.setindex!(s::SolutionODE, x, j::Int, k::Int)
    if s.n0 == 1
        @assert length(x) == 1
        @boundscheck checkbounds(s.x, j, k+1, 1)
        @inbounds setindex!(s.x, x, j, k+1, 1)
    else
        @assert length(x) == s.nd
        @boundscheck checkbounds(s.x, :, j+1, k)
        @inbounds setindex!(s.x, x, :, j+1, k)
    end
end

@inline function Base.setindex!(s::SolutionODE, x, k::Int)
    @assert length(x) == s.nd*s.n0
    @boundscheck checkbounds(s.x, :, :, k)
    @inbounds setindex!(s.x, x, :, :, k)
end



"Creates HDF5 file and initialises datasets for ODE solution object."
function createHDF5{DT,TT}(solution::SolutionODE{DT,TT}, file::AbstractString, ntime::Int=1)
    @assert ntime ≥ 1

    info("Creating HDF5 file ", file)
    # TODO Put warning if file exists.
    h5 = h5open(file, "w")

    # create dataset and copy initial conditions
    # ntime can be used to set the expected total number of timesteps
    # so that the size of the array does not need to be adapted dynamically.
    # Right now, it has to be set as dynamical size adaptation is not yet
    # working.
    x = d_create(h5, "x", datatype(DT), dataspace(solution.nd, ntime))
    x[1:solution.nd, 1] = solution[1:solution.nd, 0]

    return h5
end

"Append solution to HDF5 file."
function writeSolutionToHDF5(solution::SolutionODE, h5::HDF5.HDF5File, offset=0)
    # aquire dataset from HDF5 file
    x = h5["x"]

    # set convenience variables and compute ranges
    d  = solution.nd
    n  = solution.nt
    j1 = offset+2
    j2 = offset+1+n

    # # extend dataset if necessary
    # if size(x, 2) < j2
    #     set_dims!(x, (d, j2))
    # end

    # copy data from solution to HDF5 dataset
    x[1:d, j1:j2] = solution[1:d, 1:n]

    return nothing
end
