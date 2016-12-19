"""
`SolutionODE`: Solution of an ordinary differential equation

Contains all fields necessary to store the solution of an ODE.

### Fields

* `nd`: dimension of the dynamical variable ``q``
* `nt`: number of time steps to store
* `n0`: number of initial conditions
* `t`:  time steps
* `q`:  solution `q[nd, nt+1, n0]` with `q[:,0,:]` the initial conditions
* `ntime`: number of time steps to compute
* `nsave`: save every nsave'th time step

"""
immutable SolutionODE{dType, tType, N} <: Solution{dType, tType, N}
    nd::Int
    nt::Int
    n0::Int
    t::TimeSeries{tType}
    q::DataSeries{dType,N}
    ntime::Int
    nsave::Int

    # function SolutionODE(nd, n0, ntime, nsave, Δt)
    #     @assert dType <: Number
    #     @assert tType <: Real
    #     @assert nd > 0
    #     @assert n0 > 0
    #     @assert nsave > 0
    #     @assert ntime ≥ nsave
    #     @assert mod(ntime, nsave) == 0
    #
    #     nt = div(ntime, nsave)
    #     t = TimeSeries{tType}(nt, Δt, nsave)
    #     q = DataSeries(dType, nd, nt, n0, nsave)
    #
    #     new(nd, nt, n0, t, q, ntime, nsave)
    # end
end

function SolutionODE{DT,TT,FT}(equation::ODE{DT,TT,FT}, Δt::TT, ntime::Int, nsave::Int=1)
    N  = equation.n > 1 ? 3 : 2
    nd = equation.d
    ni = equation.n
    nt = div(ntime, nsave)

    @assert DT <: Number
    @assert TT <: Real
    @assert nd > 0
    @assert ni > 0
    @assert nsave > 0
    @assert ntime ≥ nsave
    @assert mod(ntime, nsave) == 0

    t = TimeSeries{TT}(nt, Δt, nsave)
    q = DataSeries(DT, nd, nt, ni)
    s = SolutionODE{DT,TT,N}(nd, nt, ni, t, q, ntime, nsave)
    set_initial_conditions!(s, equation)
    return s
end

function set_initial_conditions!{DT,TT}(sol::SolutionODE{DT,TT}, equ::ODE{DT,TT})
    set_initial_conditions!(sol, equ.t₀, equ.q₀)
end

function set_initial_conditions!{DT,TT}(sol::SolutionODE{DT,TT}, t₀::TT, q₀::Array{DT})
    set_data!(sol.q, q₀, 0)
    compute_timeseries!(sol.t, t₀)
end

function copy_solution!{DT,TT}(sol::SolutionODE{DT,TT}, x::Vector{DT}, n, k)
    if mod(n, sol.nsave) == 0
        set_data!(sol.q, x, div(n, sol.nsave), k)
    end
end

function reset!(s::SolutionODE)
    reset!(s.q)
    compute_timeseries!(solution.t, solution.t[end])
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
