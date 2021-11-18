
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
function Solution(equation::AbstractEquationODE, Δt, ntime::Int; kwargs...)
    SolutionODE(equation, Δt, ntime; kwargs...)
end

# Create solution for partitioned ODE.
function Solution(equation::AbstractEquationPODE, Δt, ntime::Int; kwargs...)
    SolutionPODE(equation, Δt, ntime; kwargs...)
end

# Create solution for DAE.
function Solution(equation::AbstractEquationDAE, Δt, ntime::Int; kwargs...)
    SolutionDAE(equation, Δt, ntime; kwargs...)
end

# Create solution for partitioned DAE.
function Solution(equation::AbstractEquationPDAE, Δt, ntime::Int; kwargs...)
    SolutionPDAE(equation, Δt, ntime; kwargs...)
end

# Print error for solutions of equations not implemented, yet.
function Solution(equation::Equation, Δt, ntime::Int; kwargs...)
    error("No solution found for equation ", equation)
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
# function ParallelSolution(equation::Equation, Δt, ntime::Int; kwargs...)
#     error("No parallel solution found for equation ", equation)
# end
