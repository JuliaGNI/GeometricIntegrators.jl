
using Printf

abstract type NonlinearSolver{T} end

solve!(s::NonlinearSolver) = error("solve! not implemented for $(typeof(s))")

function solve!(s::NonlinearSolver{T}, x₀::Vector{T}) where {T}
    setInitialConditions!(s, x₀)
    solve!(s)
end


struct NonlinearSolverParameters{T}
    nmin::Int   # minimum number of iterations
    nmax::Int   # maximum number of iterations
    nwarn::Int  # warn if number of iterations is larger than nwarn

    atol::T     # absolute tolerance
    rtol::T     # relative tolerance
    stol::T     # successive tolerance

    atol²::T
    rtol²::T
    stol²::T

    atol_break::T
    rtol_break::T
    stol_break::T

    function NonlinearSolverParameters{T}(nmin, nmax, nwarn, atol, rtol, stol, atol_break, rtol_break, stol_break) where {T}
        @assert nmin ≥ 0
        @assert nmax > 0
        @assert nwarn ≥ 0
        @assert atol > 0
        @assert rtol > 0
        @assert stol > 0
        @assert atol_break > 0
        @assert rtol_break > 0
        @assert stol_break > 0

        new(nmin, nmax, nwarn, atol, rtol, stol, atol^2, rtol^2, stol^2, atol_break, rtol_break, stol_break)
    end
end

function NonlinearSolverParameters(T)
    NonlinearSolverParameters{T}(get_config(:nls_nmin),
                                 get_config(:nls_nmax),
                                 get_config(:nls_nwarn),
                                 get_config(:nls_atol),
                                 get_config(:nls_rtol),
                                 get_config(:nls_stol),
                                 get_config(:nls_atol_break),
                                 get_config(:nls_rtol_break),
                                 get_config(:nls_stol_break))
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

function print_solver_status(status::NonlinearSolverStatus, params::NonlinearSolverParameters)
    if (get_config(:verbosity) == 1 && !(check_solver_converged(status, params) && status.i ≤ params.nmax)) ||
        get_config(:verbosity) > 1
        println((@sprintf "  i=%07i" status.i), ",", status)
    end
end

function check_solver_converged(status::NonlinearSolverStatus, params::NonlinearSolverParameters)
    return status.rₐ ≤ params.atol  ||
           status.rᵣ ≤ params.rtol  ||
           status.rₛ ≤ params.stol
        #    all(x -> x==0, status.Δx) ||
        #    maximum(status.Δx) / maximum(abs.(status.xₚ)) ≤ params.stol
end

function check_solver_status(status::NonlinearSolverStatus, params::NonlinearSolverParameters)
    if any(x -> isnan(x), status.xₚ)
        error("Detected NaN")
    end

    if status.rₐ > params.atol_break
        error("Absolute error of nonlinear solver ($(status.rₐ)) larger than allowed ($(params.atol_break))")
    end

    if status.rᵣ > params.rtol_break
        error("Relative error of nonlinear solver ($(status.rᵣ)) larger than allowed ($(params.rtol_break))")
    end

    if status.rₛ > params.stol_break
        error("Succesive error of nonlinear solver ($(status.rₛ)) larger than allowed ($(params.stol_break))")
    end
end


function getLinearSolver(T, n)
    linear_solver = get_config(:ls_solver)

    if linear_solver == nothing || linear_solver == :lapack
    # if linear_solver == :lapack
        linear_solver = LUSolverLAPACK{T}(BlasInt(n))
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
    status.xₚ .= x
    status.yₚ .= y
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
