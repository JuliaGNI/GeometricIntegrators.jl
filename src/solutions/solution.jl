
function Solution(problem::Union{EquationProblem, SubstepProblem}, args...; kwargs...)
    GeometricSolution(problem, args...; kwargs...)
end
    
function Solution(problem::EnsembleProblem, args...; kwargs...)
    EnsembleSolution(problem, args...; kwargs...)
end

const SolutionODE{dType, tType, dsType, probType, perType} = GeometricSolution{dType, tType, dsType, probType, perType} where {probType <: Union{ODEProblem, SODEProblem, SubstepProblem}}
const SolutionDAE{dType, tType, dsType, probType, perType} = GeometricSolution{dType, tType, dsType, probType, perType} where {probType <: DAEProblem}
const SolutionSDE{dType, tType, dsType, probType, perType} = GeometricSolution{dType, tType, dsType, probType, perType} where {probType <: SDEProblem}
const SolutionPODE{dType, tType, dsType, probType, perType} = GeometricSolution{dType, tType, dsType, probType, perType} where {probType <: Union{PODEProblem, HODEProblem, IODEProblem, LODEProblem}}
const SolutionPDAE{dType, tType, dsType, probType, perType} = GeometricSolution{dType, tType, dsType, probType, perType} where {probType <: Union{PDAEProblem, HDAEProblem, IDAEProblem, LDAEProblem}}
const SolutionPSDE{dType, tType, dsType, probType, perType} = GeometricSolution{dType, tType, dsType, probType, perType} where {probType <: PSDEProblem}


function Base.setindex!(sol::GeometricSolution, solstep::SolutionStep, n)
    sol[n] = current(solstep)
end


_periodicity(q, periodicity) = periodicity
_periodicity(q, periodicity::NullPeriodicity) = zero(q)
