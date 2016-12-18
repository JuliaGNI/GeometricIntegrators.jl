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
immutable SolutionODE{dType, tType, N} <: Solution{dType, tType, N}
    nd::Int
    nt::Int
    n0::Int
    t::Timeseries{tType}
    x::Array{dType,N}
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

        @assert N ∈ (2,3)

        if N == 2
            x = zeros(dType, nd, nt+1)
        elseif N == 3
            x = zeros(dType, nd, nt+1, n0)
        end

        new(nd, nt, n0, t, x, ntime, nsave)
    end
end

function SolutionODE{DT,TT,FT}(equation::ODE{DT,TT,FT}, Δt::TT, ntime::Int, nsave::Int=1)
    N = equation.n > 1 ? 3 : 2
    s = SolutionODE{DT,TT,N}(equation.d, equation.n, ntime, nsave, Δt)
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

function reset{DT,TT}(s::SolutionODE{DT,TT,2})
    for i in 1:size(solution,1)
        solution[i, 0] = solution[i, end]
    end
    solution.t[0] = solution.t[end]
    compute_timeseries!(solution.t)
end

function reset{DT,TT}(s::SolutionODE{DT,TT,3})
    for k in 1:size(solution,3)
        for i in 1:size(solution,1)
            solution[i, 0, k] = solution[i, end, k]
        end
    end
    solution.t[0] = solution.t[end]
    compute_timeseries!(solution.t)
end

Base.indices{DT,TT}(s::SolutionODE{DT,TT,2}) = (1:s.nd, 0:s.nt)
Base.indices{DT,TT}(s::SolutionODE{DT,TT,3}) = (1:s.nd, 0:s.nt, 1:s.n0)
Base.strides(s::SolutionODE) = strides(s.x)

@inline function Base.getindex{DT,TT}(s::SolutionODE{DT,TT,2}, i::Int, j::Int)
    @boundscheck checkbounds(s.x, i, j+1)
    @inbounds r = getindex(s.x, i, j+1)
    return r
end

@inline function Base.getindex{DT,TT}(s::SolutionODE{DT,TT,2}, j::Int)
    @boundscheck checkbounds(s.x, :, j+1)
    @inbounds r = getindex(s.x, :, j+1)
    return r
end

@inline function Base.getindex{DT,TT}(s::SolutionODE{DT,TT,3}, i::Int, j::Int, k::Int)
    @boundscheck checkbounds(s.x, i, j+1, k)
    @inbounds r = getindex(s.x, i, j+1, k)
    return r
end

@inline function Base.getindex{DT,TT}(s::SolutionODE{DT,TT,3}, j::Int, k::Int)
    @boundscheck checkbounds(s.x, :, j+1, k)
    @inbounds r = getindex(s.x, :, j+1, k)
    return r
end

@inline function Base.getindex{DT,TT}(s::SolutionODE{DT,TT,3}, k::Int)
    @boundscheck checkbounds(s.x, :, :, k)
    @inbounds r = getindex(s.x, :, :, k)
    return r
end

@inline function Base.setindex!{DT,TT}(s::SolutionODE{DT,TT,2}, x, i::Int, j::Int)
    @assert length(x) == 1
    @boundscheck checkbounds(s.x, i, j+1)
    @inbounds setindex!(s.x, x, i, j+1)
end

@inline function Base.setindex!{DT,TT}(s::SolutionODE{DT,TT,2}, x, j::Int)
    @assert ndims(x) == 1
    @assert length(x) == size(s.x, 1)
    @boundscheck checkbounds(s.x, :, j)
    @inbounds setindex!(s.x, x, :, j)
end

@inline function Base.setindex!{DT,TT}(s::SolutionODE{DT,TT,3}, x, i::Int, j::Int, k::Int)
    @assert length(x) == 1
    @boundscheck checkbounds(s.x, i, j+1, k)
    @inbounds setindex!(s.x, x, i, j+1, k)
end

@inline function Base.setindex!{DT,TT}(s::SolutionODE{DT,TT,3}, x, j::Int, k::Int)
    @assert ndims(x) == 1
    @assert length(x) == size(s.x, 1)
    @boundscheck checkbounds(s.x, :, j+1, k)
    @inbounds setindex!(s.x, x, :, j+1, k)
end

@inline function Base.setindex!{DT,TT}(s::SolutionODE{DT,TT,3}, x, k::Int)
    @assert ndims(x) == 2
    @assert size(x,1) == size(s.x, 1)
    @assert size(x,2) == size(s.x, 2)
    @boundscheck checkbounds(s.x, :, :, k)
    @inbounds setindex!(s.x, x, :, :, k)
end

# TODO Fix HDF5 functions.

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
