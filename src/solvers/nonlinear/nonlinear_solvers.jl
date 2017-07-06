
abstract type NonlinearSolver{T} end

solve!(s::NonlinearSolver) = error("solve! not implemented for $(typeof(s))")

function solve!(s::NonlinearSolver{T}, x₀::Vector{T}) where {T}
    setInitialConditions!(s, x₀)
    solve!(s)
end


# default parameters for nonlinear solvers
const DEFAULT_nwarn=100

abstract type NonlinearFunctionParameters{T} end

function function_stages!(y::Vector{DT}, b::Vector{TT}, params::NonlinearFunctionParameters{DT}) where {DT,TT}
    error("No function_stages implemented for this integrator.")
end


struct NonlinearSolverParameters{T}
    nmin::Int   # minimum number of iterations
    nmax::Int   # maximum number of iterations

    atol::T     # absolute tolerance
    rtol::T     # relative tolerance
    stol::T     # successive tolerance

    atol²::T
    rtol²::T
    stol²::T

    function NonlinearSolverParameters{T}(nmin, nmax, atol, rtol, stol) where {T}
        @assert nmin ≥ 0
        @assert nmax > 0
        @assert atol > 0
        @assert rtol > 0
        @assert stol > 0

        new(nmin, nmax, atol, rtol, stol, atol^2, rtol^2, stol^2)
    end
end

function NonlinearSolverParameters(T)
    NonlinearSolverParameters{T}(get_config(:nls_nmin),
                                 get_config(:nls_nmax),
                                 get_config(:nls_atol),
                                 get_config(:nls_rtol),
                                 get_config(:nls_stol))
end


mutable struct NonlinearSolverStatus{T}
    i::Int         # iteration number
    rₐ::T          # residual (absolute)
    rᵣ::T          # residual (relative)
    rₛ::T          # residual (successive)
    r₀::Vector{T}  # initial residual (absolute)
    y₀::Vector{T}  # initial solution
    yₚ::Vector{T}  # previous solution

    NonlinearSolverStatus{T}(n) where {T} = new(0, 0, 0, 0, zeros(T,n), zeros(T,n), zeros(T,n))
end

Base.show(io::IO, status::NonlinearSolverStatus) = print(io,
                        (@sprintf "%4i"    status.i),  ", ", (@sprintf "%14.8e" status.rₐ), ", ",
                        (@sprintf "%14.8e" status.rᵣ), ", ", (@sprintf "%14.8e" status.rₛ))

function solverConverged(status::NonlinearSolverStatus, params::NonlinearSolverParameters)
    return status.rₐ ≤ params.atol² ||
           status.rᵣ ≤ params.rtol²# ||
        #    status.rₛ ≤ params.stol²
end

function solverStatusOK(status::NonlinearSolverStatus, params::NonlinearSolverParameters)
    return solverConverged(status, params) && status.i  ≤ params.nmax
end

function getLinearSolver(T, n)
    linear_solver = get_config(:ls_solver)

    if linear_solver == nothing || linear_solver == :lapack
    # if linear_solver == :lapack
        linear_solver = LUSolverLAPACK{T}(n)
    elseif linear_solver == :julia
        linear_solver = LUSolver{T}(n)
    # elseif linear_solver == nothing
    #     linear_solver = IterativeSolver{T}(n)
    else
        @assert typeof(linear_solver) <: LinearSolver{T}
        @assert n == linear_solver.n
    end
    return linear_solver
end


function residual_initial!(status::NonlinearSolverStatus{T}, y::Vector{T}) where {T}
    @assert length(y) == length(status.r₀)
    @inbounds for i in 1:length(status.r₀)
        status.r₀[i] = y[i]^2
    end
    status.rₐ = maximum(status.r₀)
    status.rᵣ = 1
    simd_copy!(y, status.y₀)
    simd_copy!(y, status.yₚ)
end

function residual!(status::NonlinearSolverStatus{T}, y::Vector{T}) where {T}
    residual_absolute!(status, y)
    residual_relative!(status, y)
    residual_successive!(status, y)
    simd_copy!(y, status.yₚ)
end

function residual_absolute!(status::NonlinearSolverStatus{T}, y::Vector{T}) where {T}
    @assert length(y) == length(status.r₀)
    status.rₐ = 0
    @inbounds for yᵢ in y
        status.rₐ = max(status.rₐ, yᵢ^2)
    end
end

function residual_relative!(status::NonlinearSolverStatus{T}, y::Vector{T}) where {T}
    @assert length(y) == length(status.r₀)
    status.rᵣ = 0
    @inbounds for i in 1:length(y)
        status.rᵣ = max(status.rᵣ, y[i]^2 / status.r₀[i])
    end
end

function residual_successive!(status::NonlinearSolverStatus{T}, y::Vector{T}) where {T}
    @assert length(y) == length(status.yₚ)
    status.rₛ = 0
    @inbounds for i in 1:length(y)
        status.rₛ = max(status.rₛ, ( (status.yₚ[i] - y[i]) / status.y₀[i] )^2)
    end
end
