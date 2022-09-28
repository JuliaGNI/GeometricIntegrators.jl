
const Solution = GeometricSolution

const SolutionODE{dType, tType, dsType, probType, perType} = GeometricSolution{dType, tType, dsType, probType, perType} where {probType <: AbstractProblemODE}
const SolutionDAE{dType, tType, dsType, probType, perType} = GeometricSolution{dType, tType, dsType, probType, perType} where {probType <: AbstractProblemDAE}
const SolutionSDE{dType, tType, dsType, probType, perType} = GeometricSolution{dType, tType, dsType, probType, perType} where {probType <: AbstractProblemSDE}
const SolutionPODE{dType, tType, dsType, probType, perType} = GeometricSolution{dType, tType, dsType, probType, perType} where {probType <: AbstractProblemPODE}
const SolutionPDAE{dType, tType, dsType, probType, perType} = GeometricSolution{dType, tType, dsType, probType, perType} where {probType <: AbstractProblemPDAE}
const SolutionPSDE{dType, tType, dsType, probType, perType} = GeometricSolution{dType, tType, dsType, probType, perType} where {probType <: AbstractProblemPSDE}


function Base.setindex!(sol::Solution, solstep::SolutionStep, n)
    sol[n] = current(solstep)
end


_periodicity(q, periodicity) = periodicity
_periodicity(q, periodicity::NullPeriodicity) = zero(q)
