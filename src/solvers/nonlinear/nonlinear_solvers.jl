
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

    x₀::Vector{T}  # initial solution
    xₚ::Vector{T}  # previous solution

    y₀::Vector{T}  # initial function
    yₚ::Vector{T}  # previous function

    NonlinearSolverStatus{T}(n) where {T} = new(0, 0, 0, 0, zeros(T,n),
                                                zeros(T,n), zeros(T,n),
                                                zeros(T,n), zeros(T,n),)
end

Base.show(io::IO, status::NonlinearSolverStatus) = print(io,
                        (@sprintf "    n=%4i" status.i),  ",   ", (@sprintf "rₐ=%14.8e" status.rₐ), ",   ",
                        (@sprintf "rᵣ=%14.8e" status.rᵣ), ",   ", (@sprintf "rₛ=%14.8e" status.rₛ))

function printSolverStatus(status::NonlinearSolverStatus, params::NonlinearSolverParameters, n::Int)
    if (get_config(:verbosity) == 1 && !solverStatusOK(status, params)) ||
        get_config(:verbosity) > 1
        println((@sprintf "  i=%07i" n), ",", status)
    end
end

function printSolverStatus(status::NonlinearSolverStatus, params::NonlinearSolverParameters)
    if (get_config(:verbosity) == 1 && !solverStatusOK(status, params)) ||
        get_config(:verbosity) > 1
        println(status)
    end
end

function solverConverged(status::NonlinearSolverStatus, params::NonlinearSolverParameters)
    return status.rₐ ≤ params.atol  ||
           status.rᵣ ≤ params.rtol  ||
           status.rₛ ≤ params.stol
        #    all(x -> x==0, status.Δx) ||
        #    maximum(status.Δx) / maximum(abs.(status.xₚ)) ≤ params.stol
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


function residual_initial!(status::NonlinearSolverStatus{T}, x::Vector{T}, y::Vector{T}) where {T}
    @assert length(x) == length(y) == length(status.r₀)

    status.r₀ .= y.^2
    status.rₐ  = maximum(status.r₀)
    status.rᵣ  = 1
    status.x₀ .= x
    status.xₚ .= x
    status.y₀ .= y
    status.yₚ .= y
end

function residual!(status::NonlinearSolverStatus{T}, δx::Vector{T}, x::Vector{T}, y::Vector{T}) where {T}
    residual_absolute!(status, y)
    residual_relative!(status, y)
    residual_successive!(status, δx, x)
    simd_copy!(x, status.xₚ)
    simd_copy!(y, status.yₚ)
end

function residual_absolute!(status::NonlinearSolverStatus{T}, y::Vector{T}) where {T}
    @assert length(y) == length(status.y₀)
    local rₐ²::T = 0
    @inbounds for yᵢ in y
        rₐ² = max(rₐ², yᵢ^2)
    end
    status.rₐ = sqrt(rₐ²)
end

function residual_relative!(status::NonlinearSolverStatus{T}, y::Vector{T}) where {T}
    @assert length(y) == length(status.y₀)
    status.rᵣ = norm(y) / norm(status.y₀)
end

function residual_successive!(status::NonlinearSolverStatus{T}, δx::Vector{T}, x::Vector{T}) where {T}
    @assert length(δx) == length(x)
    status.rₛ = norm(δx) / norm(x)
end
