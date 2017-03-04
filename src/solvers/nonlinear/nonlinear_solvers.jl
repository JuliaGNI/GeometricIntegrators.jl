
abstract NonlinearSolver{T}

solve!(s::NonlinearSolver) = error("solve! not implemented for $(typeof(s))")

function solve!{T}(s::NonlinearSolver{T}, x₀::Vector{T})
    setInitialConditions!(s, x₀)
    solve!(s)
end


# default parameters for nonlinear solvers
const DEFAULT_atol=2eps()
const DEFAULT_rtol=2eps()
const DEFAULT_stol=2eps()
const DEFAULT_nmax=1000
const DEFAULT_nwarn=100
const DEFAULT_ϵ=4sqrt(eps())

abstract NonlinearFunctionParameters{T}

function function_stages!{DT,TT}(y::Vector{DT}, b::Vector{TT}, params::NonlinearFunctionParameters{DT})
    error("No function_stages implemented for this integrator.")
end


immutable NonlinearSolverParameters{T}
    nmax::Int   # maximum number of iterations

    atol::T     # absolute tolerance
    rtol::T     # relative tolerance
    stol::T     # successive tolerance

    atol²::T
    rtol²::T
    stol²::T

    function NonlinearSolverParameters(nmax, atol, rtol, stol)
        @assert nmax > 0
        @assert atol > 0
        @assert rtol > 0
        @assert stol > 0

        new(nmax, atol, rtol, stol, atol^2, rtol^2, stol^2)
    end
end


type NonlinearSolverStatus{T}
    i::Int         # iteration number
    rₐ::T          # residual (absolute)
    rᵣ::T          # residual (relative)
    rₛ::T          # residual (successive)
    r₀::Vector{T}  # initial residual (absolute)
    y₀::Vector{T}  # initial solution
    yₚ::Vector{T}  # previous solution

    NonlinearSolverStatus(n) = new(0, 0, 0, 0, zeros(T,n), zeros(T,n), zeros(T,n))
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

function getLinearSolver(T, n, linear_solver)
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


function residual_initial!{T}(status::NonlinearSolverStatus{T}, y::Vector{T})
    @assert length(y) == length(status.r₀)
    @inbounds for i in 1:length(status.r₀)
        status.r₀[i] = y[i]^2
    end
    status.rₐ = maximum(status.r₀)
    status.rᵣ = 1
    simd_copy!(y, status.y₀)
    simd_copy!(y, status.yₚ)
end

function residual!{T}(status::NonlinearSolverStatus{T}, y::Vector{T})
    residual_absolute!(status, y)
    residual_relative!(status, y)
    residual_successive!(status, y)
    simd_copy!(y, status.yₚ)
end

function residual_absolute!{T}(status::NonlinearSolverStatus{T}, y::Vector{T})
    @assert length(y) == length(status.r₀)
    status.rₐ = 0
    @inbounds for yᵢ in y
        status.rₐ = max(status.rₐ, yᵢ^2)
    end
end

function residual_relative!{T}(status::NonlinearSolverStatus{T}, y::Vector{T})
    @assert length(y) == length(status.r₀)
    status.rᵣ = 0
    @inbounds for i in 1:length(y)
        status.rᵣ = max(status.rᵣ, y[i]^2 / status.r₀[i])
    end
end

function residual_successive!{T}(status::NonlinearSolverStatus{T}, y::Vector{T})
    @assert length(y) == length(status.yₚ)
    status.rₛ = 0
    @inbounds for i in 1:length(y)
        status.rₛ = max(status.rₛ, ( (status.yₚ[i] - y[i]) / status.y₀[i] )^2)
    end
end
