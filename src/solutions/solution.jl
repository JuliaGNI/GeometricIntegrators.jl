
abstract type DeterministicSolution{dType, tType, N} <: AbstractSolution{dType, tType, N} end


function parameters(solution::DeterministicSolution)
    (
        ntime = ntime(solution),
        nsave = nsave(solution),
        nsamples = nsamples(solution),
    )
end


"""
```julia
Solution(problem; kwargs...)
```

Create the appropriate `Solution` for the given `problem` type for a
simulation with `ntime` time steps of step size `Δt`.

"""
function Solution end


# Create solution for ODE.
function Solution(problem::AbstractProblemODE; kwargs...)
    SolutionODE(problem; kwargs...)
end

# Create solution for partitioned ODE.
function Solution(problem::AbstractProblemPODE; kwargs...)
    SolutionPODE(problem; kwargs...)
end

# Create solution for DAE.
function Solution(problem::AbstractProblemDAE; kwargs...)
    SolutionDAE(problem; kwargs...)
end

# Create solution for partitioned DAE.
function Solution(problem::AbstractProblemPDAE; kwargs...)
    SolutionPDAE(problem; kwargs...)
end

# Print error for solutions of problems not implemented, yet.
function Solution(problem::GeometricProblem; kwargs...)
    error("No solution found for problem type ", typeof(problem))
end


# """
# ```julia
# ParallelSolution(problem; kwargs...)
# ```

# Create the appropriate `ParallelSolution` for the given `problem` type for a
# simulation with `ntime` time steps of step size `Δt`.

# """
# function ParallelSolution end

# # Create parallel solution for ODE.
# function ParallelSolution(problem::AbstractProblemODE; kwargs...)
#     PSolutionODE(problem; kwargs...)
# end

# # Create parallel solution for partitioned ODE.
# function ParallelSolution(problem::AbstractProblemPODE; kwargs...)
#     PSolutionPODE(problem; kwargs...)
# end

# # Create parallel solution for DAE.
# function ParallelSolution(problem::AbstractProblemDAE; kwargs...)
#     PSolutionDAE(problem; kwargs...)
# end

# # Create parallel solution for partitioned DAE.
# function ParallelSolution(problem::AbstractProblemPDAE; kwargs...)
#     PSolutionPDAE(problem; kwargs...)
# end

# # Print error for parallel solutions of problems not implemented, yet.
# function ParallelSolution(problem::GeometricProblem; kwargs...)
#     error("No parallel solution found for problem ", problem)
# end


_periodicity(q, periodicity) = periodicity
_periodicity(q, periodicity::NullPeriodicity) = zero(q)
