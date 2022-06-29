
abstract type DeterministicSolution{dType, tType, N} <: Solution{dType, tType, N} end


function parameters(solution::DeterministicSolution)
    (
        ntime = ntime(solution),
        nsave = nsave(solution),
        nsamples = nsamples(solution),
    )
end


"""
```julia
Solution(equation, Δt, ntime; kwargs...)
```

Create the appropriate `Solution` for the given `equation` type for a
simulation with `ntime` time steps of step size `Δt`.

"""
function Solution end


# Create solution for ODE.
function Solution(problem::AbstractProblemODE, ntime::Int; kwargs...)
    SolutionODE(problem, problem.tstep, ntime; kwargs...)
end

# Create solution for partitioned ODE.
function Solution(problem::AbstractProblemPODE, ntime::Int; kwargs...)
    SolutionPODE(problem, problem.tstep, ntime; kwargs...)
end

# Create solution for DAE.
function Solution(problem::AbstractProblemDAE, ntime::Int; kwargs...)
    SolutionDAE(problem, problem.tstep, ntime; kwargs...)
end

# Create solution for partitioned DAE.
function Solution(problem::AbstractProblemPDAE, ntime::Int; kwargs...)
    SolutionPDAE(problem, problem.tstep, ntime; kwargs...)
end

# Print error for solutions of equations not implemented, yet.
function Solution(problem::GeometricProblem, ntime::Int; kwargs...)
    error("No solution found for problem type ", typeof(problem))
end


# """
# ```julia
# ParallelSolution(equation, Δt, ntime; kwargs...)
# ```

# Create the appropriate `ParallelSolution` for the given `equation` type for a
# simulation with `ntime` time steps of step size `Δt`.

# """
# function ParallelSolution end

# # Create parallel solution for ODE.
# function ParallelSolution(equation::AbstractEquationODE, Δt, ntime::Int; kwargs...)
#     PSolutionODE(equation, Δt, ntime; kwargs...)
# end

# # Create parallel solution for partitioned ODE.
# function ParallelSolution(equation::AbstractEquationPODE, Δt, ntime::Int; kwargs...)
#     PSolutionPODE(equation, Δt, ntime; kwargs...)
# end

# # Create parallel solution for DAE.
# function ParallelSolution(equation::AbstractEquationDAE, Δt, ntime::Int; kwargs...)
#     PSolutionDAE(equation, Δt, ntime; kwargs...)
# end

# # Create parallel solution for partitioned DAE.
# function ParallelSolution(equation::AbstractEquationPDAE, Δt, ntime::Int; kwargs...)
#     PSolutionPDAE(equation, Δt, ntime; kwargs...)
# end

# # Print error for parallel solutions of equations not implemented, yet.
# function ParallelSolution(equation::GeometricEquation, Δt, ntime::Int; kwargs...)
#     error("No parallel solution found for equation ", equation)
# end


_periodicity(q, periodicity) = periodicity
_periodicity(q, periodicity::NullPeriodicity) = zero(q)
