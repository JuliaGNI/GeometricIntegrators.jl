


const Solution = GeometricSolution

const SolutionODE{dType, tType, dsType, probType, perType} = GeometricSolution{dType, tType, dsType, probType, perType} where {probType <: AbstractProblemODE}
const SolutionDAE{dType, tType, dsType, probType, perType} = GeometricSolution{dType, tType, dsType, probType, perType} where {probType <: AbstractProblemDAE}
const SolutionSDE{dType, tType, dsType, probType, perType} = GeometricSolution{dType, tType, dsType, probType, perType} where {probType <: AbstractProblemSDE}
const SolutionPODE{dType, tType, dsType, probType, perType} = GeometricSolution{dType, tType, dsType, probType, perType} where {probType <: AbstractProblemPODE}
const SolutionPDAE{dType, tType, dsType, probType, perType} = GeometricSolution{dType, tType, dsType, probType, perType} where {probType <: AbstractProblemPDAE}
const SolutionPSDE{dType, tType, dsType, probType, perType} = GeometricSolution{dType, tType, dsType, probType, perType} where {probType <: AbstractProblemPSDE}


function Base.setindex!(sol::Solution, asol::AtomicSolution, n)
    sol[n] = current(asol)
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
