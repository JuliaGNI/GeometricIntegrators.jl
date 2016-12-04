
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
const DEFAULT_nmax=100
const DEFAULT_ϵ=sqrt(eps())


abstract NonlinearFunctionParameters{T}

function function_stages!{T}(y::Vector{T}, b::Vector{T}, params::NonlinearFunctionParameters{T})
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
    i::Int      # iteration number
    r₀::T       # initial residual (absolute)
    rₐ::T       # residual (absolute)
    rᵣ::T       # residual (relative)
    rₛ::T       # residual (successive)

    NonlinearSolverStatus() = new(0, 0., 0., 0., 0.)
end

Base.show(io::IO, status::NonlinearSolverStatus) = print(io,
                        status.i, ", ", status.rₐ,", ",
                        status.rᵣ,", ", status.rₛ)

function solverStatusOK(status::NonlinearSolverStatus, params::NonlinearSolverParameters)
    return status.rₐ < params.atol² &&
           status.rᵣ < params.rtol  &&
           status.rₛ < params.stol² &&
           status.i  ≤ params.nmax
end

function getLinearSolver(T, n, linear_solver)
    if linear_solver == nothing
        linear_solver = LUSolverLAPACK(zeros(T, n, n), zeros(T, n))
    else
        @assert typeof(linear_solver) <: LinearSolver{T}
        @assert n == linear_solver.n
    end
    return linear_solver
end


function residual_absolute{T}(x::Vector{T})
    local r::T = 0.
    for xᵢ in x
        r = max(r, xᵢ*xᵢ)
    end
    r
end

function residual_relative{T}(δx::Vector{T}, x::Vector{T})
    @assert length(x) == length(δx)
    local r::T = 0.
    for i in 1:length(x)
        r = max(r, abs(δx[i] / x[i]))
    end
    r
end
