
using HDF5

abstract Solution{T,N} <: DenseArray{T,N}

"Create solution for ODE."
function Solution(equation::ODE, ntime::Int, nsave::Int=1)
    SolutionODE(equation, ntime, nsave)
end

"Create solution for partitioned ODE."
function Solution(equation::PODE, ntime::Int, nsave::Int=1)
    SolutionPODE(equation, ntime, nsave)
end

"Create solution for special partitioned ODE."
function Solution(equation::SPODE, ntime::Int, nsave::Int=1)
    SolutionPODE(equation, ntime, nsave)
end

"Create solution for DAE."
function Solution(equation::DAE, ntime::Int, nsave::Int=1)
    SolutionDAE(equation, ntime, nsave)
end

"Create solution for partitioned DAE."
function Solution(equation::PDAE, ntime::Int, nsave::Int=1)
    SolutionPDAE(equation, ntime, nsave)
end

"Print error for solutions of equations not implemented, yet."
function Solution(equation::Equation, ntime::Int, nsave::Int=1)
    error("No solution found for equation ", equation)
end

# "createHDF5: Creates or opens HDF5 file."
# function createHDF5(file::AbstractString, overwrite=true)
#     if overwrite
#         flag = "w"
#     else
#         flag = "r+"
#     end
#
#     h5open(file, flag)
# end

"Creates HDF5 file, writes solution to file, and closes file."
function writeSolutionToHDF5(solution::Solution, file::AbstractString)
    h5 = createHDF5(solution, file, solution.n+1)
    writeSolutionToHDF5(solution, h5)
    close(h5)
end


Base.eltype{T,N}(s::Solution{T,N}) = T
Base.ndims{T,N}(s::Solution{T,N}) = N
Base.size(s::Solution) = size(s.x)
# Base.length(s::Solution) = s.d * s.n
# Base.length(s::Solution) = length(s.x)
# Base.endof(s::Solution) = length(s)
Base.indices(s::Solution, d) = indices(s)[d]
Base.stride(s::Solution, d) = strides(s)[d]

# TODO Implement similar() and convert() to/from array functions.


"Solution of an ordinary differential equation."
immutable SolutionODE{T} <: Solution{T,2}
    d::Int
    n::Int
    t::Array{T,1}
    x::Array{T,2}
    ntime::Int
    nsave::Int

    function SolutionODE(d::Int, ntime::Int, nsave::Int)
        @assert T <: Real
        @assert d > 0
        @assert nsave > 0
        @assert ntime ≥ nsave
        @assert mod(ntime, nsave) == 0

        n = div(ntime, nsave)
        t = zeros(T, n+1)
        x = zeros(T, d, n+1)
        new(d, n, t, x, ntime, nsave)
    end
end

function SolutionODE{T}(equation::ODE{T}, ntime::Int, nsave::Int=1)
    s = SolutionODE{T}(equation.d, ntime, nsave)
    set_initial_conditions!(s, equation)
    return s
end

function set_initial_conditions!(solution::SolutionODE, equation::ODE)
    simd_copy_yx_first!(equation.q₀, solution, 0)
    solution.t[1] = equation.t₀
end

function reset(s::SolutionODE)
    solution[1:solution.d, 0] = solution[1:solution.d, solution.n]
    solution.t[0] = equation.t[solution.n+1]
end

Base.indices(s::SolutionODE) = (1:s.d, 0:s.n)
Base.strides(s::SolutionODE) = (1, s.d)

@inline function Base.getindex(s::SolutionODE, i::Int, j::Int)
    @boundscheck checkbounds(s.x, i, j+1)
    @inbounds r = getindex(s.x, i, j+1)
    return r
end

@inline function Base.getindex(s::SolutionODE, j::Int)
    @boundscheck checkbounds(s.x, 1:s.d, j+1)
    @inbounds r = getindex(s.x, 1:s.d, j+1)
    return r
end

@inline function Base.setindex!(s::SolutionODE, x, i::Int, j::Int)
    @boundscheck checkbounds(s.x, i, j+1)
    @inbounds setindex!(s.x, x, i, j+1)
end

@inline function Base.setindex!(s::SolutionODE, x, j::Int)
    @assert length(x) == s.d
    @boundscheck checkbounds(s.x, 1:s.d, j+1)
    @inbounds setindex!(s.x, x, 1:s.d, j+1)
end



"Creates HDF5 file and initialises datasets for ODE solution object."
function createHDF5{T}(solution::SolutionODE{T}, file::AbstractString, ntime::Int=1)
    @assert ntime ≥ 1

    info("Creating HDF5 file ", file)
    # TODO Put warning if file exists.
    h5 = h5open(file, "w")

    # create dataset and copy initial conditions
    # ntime can be used to set the expected total number of timesteps
    # so that the size of the array does not need to be adapted dynamically.
    # Right now, it has to be set as dynamical size adaptation is not yet
    # working.
    x = d_create(h5, "x", datatype(T), dataspace(solution.d, ntime))
    x[1:solution.d, 1] = solution[1:solution.d, 0]

    return h5
end

"Append solution to HDF5 file."
function writeSolutionToHDF5(solution::SolutionODE, h5::HDF5.HDF5File, offset=0)
    # aquire dataset from HDF5 file
    x = h5["x"]

    # set convenience variables and compute ranges
    d  = solution.d
    n  = solution.n
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


"Solution of a partitioned ordinary differential equation."
immutable SolutionPODE{T} <: Solution{T,3}
    d::Int
    n::Int
    t::Array{T,1}
    x::Array{T,3}
    q::AbstractArray{T,2}
    p::AbstractArray{T,2}
    ntime::Int
    nsave::Int

    function SolutionPODE(d, ntime, nsave)
        @assert T <: Real
        @assert d > 0
        @assert nsave > 0
        @assert ntime ≥ nsave
        @assert mod(ntime, nsave) == 0

        n = div(ntime, nsave)
        t = zeros(T, n+1)
        x = zeros(T, d, 2, n+1)
        q = view(x, :, 1, 1:n+1)
        p = view(x, :, 2, 1:n+1)
        new(d, n, t, x, q, p, ntime, nsave)
    end
end

function SolutionPODE(equation::PODE, ntime::Int, nsave::Int=1)
    T1 = eltype(equation.q₀)
    T2 = eltype(equation.p₀)
    @assert T1 == T2
    s = SolutionPODE{T1}(equation.d, ntime, nsave)
    set_initial_conditions!(s, equation)
    return s
end

function SolutionPODE(equation::SPODE, ntime::Int, nsave::Int=1)
    T1 = eltype(equation.q₀)
    T2 = eltype(equation.p₀)
    @assert T1 == T2
    s = SolutionPODE{T1}(equation.d, ntime, nsave)
    set_initial_conditions!(s, equation)
    return s
end

function set_initial_conditions!(solution::SolutionPODE, equation::PODE)
    solution[1:solution.d, 1, 0] = equation.q₀
    solution[1:solution.d, 2, 0] = equation.p₀
end

function set_initial_conditions!(solution::SolutionPODE, equation::SPODE)
    solution[1:solution.d, 1, 0] = equation.q₀
    solution[1:solution.d, 2, 0] = equation.p₀
end

function reset!(s::SolutionPODE)
    solution[1:solution.d, 1:2, 0] = solution[1:solution.d, 1:2, solution.n]
end

Base.indices(s::SolutionPODE) = (1:s.d, 1:2, 0:s.n)
Base.strides(s::SolutionPODE) = (1, s.d, 2*s.d)

@inline function Base.getindex(s::SolutionPODE, i::Int, j::Int, k::Int)
    @boundscheck checkbounds(s.x, i, j, k+1)
    @inbounds r = getindex(s.x, i, j, k+1)
    return r
end

@inline function Base.getindex(s::SolutionPODE, j::Int, k::Int)
    @boundscheck checkbounds(s.x, 1:s.d, j, k+1)
    @inbounds r = getindex(s.x, 1:s.d, j, k+1)
    return r
end

@inline function Base.setindex!(s::SolutionPODE, x, i::Int, j::Int, k::Int)
    @boundscheck checkbounds(s.x, i, j, k+1)
    @inbounds setindex!(s.x, x, i, j, k+1)
end

@inline function Base.setindex!(s::SolutionPODE, x, j::Int, k::Int)
    @assert length(x) == s.d
    @boundscheck checkbounds(s.x, 1:s.d, j, k+1)
    @inbounds setindex!(s.x, x, 1:s.d, j, k+1)
end


"Solution of a differential algebraic equation."
immutable SolutionDAE{T} <: Solution{T,3}
    d::Int
    m::Int
    n::Int
    t::Array{T,1}
    x::Array{T,2}
    λ::Array{T,2}
    ntime::Int
    nsave::Int

    function SolutionDAE(d, m, ntime, nsave)
        @assert T <: Real
        @assert d > 0
        @assert m > 0
        @assert nsave > 0
        @assert ntime ≥ nsave
        @assert mod(ntime, nsave) == 0

        n = div(ntime, nsave)
        t = zeros(T, n+1)
        x = zeros(T, d, n+1)
        λ = zeros(T, m, n+1)
        new(d, m, n, t, x, λ, ntime, nsave)
    end
end

function SolutionDAE(equation::DAE, ntime::Int, nsave::Int=1)
    T1 = eltype(equation.q₀)
    T2 = eltype(equation.λ₀)
    @assert T1 == T2
    SolutionDAE{T1}(equation.m, equation.n, ntime, nsave)
end

function Base.getindex(s::SolutionDAE, i::Int)
    # TODO
end

function reset(s::SolutionDAE)
    # TODO
end


"Solution of a partitioned differential algebraic equation."
immutable SolutionPDAE{T} <: Solution{T,3}
    d::Int
    m::Int
    n::Int
    t::Array{T,1}
    q::Array{T,2}
    p::Array{T,2}
    λ::Array{T,2}
    ntime::Int
    nsave::Int

    function SolutionPDAE(d, m, ntime, nsave)
        @assert T <: Real
        @assert d > 0
        @assert m > 0
        @assert nsave > 0
        @assert ntime ≥ nsave
        @assert mod(ntime, nsave) == 0

        n = div(ntime, nsave)
        t = zeros(T, n+1)
        q = zeros(T, d, n+1)
        p = zeros(T, d, n+1)
        λ = zeros(T, m, n+1)
        new(d, m, n, t, q, p, λ, ntime, nsave)
    end
end

function SolutionPDAE(equation::PDAE, ntime::Int, nsave::Int=1)
    T1 = eltype(equation.q₀)
    T2 = eltype(equation.p₀)
    T3 = eltype(equation.λ₀)
    @assert T1 == T2 == T3
    SolutionPDAE{T1}(equation.m, equation.n, ntime, nsave)
end

function Base.getindex(s::SolutionPDAE, i::Int)
    # TODO
end

function reset(s::SolutionPDAE)
    # TODO
end
