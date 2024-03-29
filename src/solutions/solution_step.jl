
default_extrapolation() = MidpointExtrapolation(5)

"""
Abstract atomic or single-step solution.

## Constructors:

```julia
SolutionStep(equation)
SolutionStep(solution)
SolutionStep(solution, integrator)
```

Automatically construct the appropriate atomic solution based on the 
given `equation` or `solution` type. If an `integrator` is provided as,
the `internal` field of the atomic solution is constructed according to
the internal state of the integrator as obtained from the function
`internal_variables`.

"""
abstract type SolutionStep{dType <: Number, tType <: Real, aType <: AbstractArray{dType}} end

initialize!(::SolutionStep, ::AbstractProblem) = nothing

GeometricBase.datatype(::SolutionStep{dType, tType, aType}) where {dType, tType, aType} = dType
GeometricBase.timetype(::SolutionStep{dType, tType, aType}) where {dType, tType, aType} = tType
GeometricBase.arrtype(::SolutionStep{dType, tType, aType}) where {dType, tType, aType} = aType

eachhistory(sol::SolutionStep) = 1:nhistory(sol)
backwardhistory(sol::SolutionStep) = nhistory(sol):-1:1


"""
Copy the initial conditions of a `EquationProblem` to the current state of an atomic solution.
"""
function Base.copy!(solstep::SolutionStep, equ::EquationProblem)
    copy!(solstep, initial_conditions(equ))
end

"""
Returns a NamedTuple with the solution of the current time step.
"""
function current end

"""
Returns a NamedTuple with the solution of the previous time step.
"""
function previous end

function cut_periodic_solution!(solstep::SolutionStep, periodicity)
    if periodicity.q != NullPeriodicity()
        @assert axes(solstep.q) == axes(periodicity.q)
        for k in eachindex(solstep.q, periodicity.q)
            if periodicity.q[k] ≠ 0
                while solstep.q[k] < 0
                    solstep.q[k], solstep.q̃[k] = compensated_summation(+periodicity.q[k], solstep.q[k], solstep.q̃[k])
                    for i in eachhistory(solstep)
                        history(solstep).q[i][k] += periodicity.q[k]
                    end
                end
                while solstep.q[k] ≥ periodicity.q[k]
                    solstep.q[k], solstep.q̃[k] = compensated_summation(-periodicity.q[k], solstep.q[k], solstep.q̃[k])
                    for i in eachhistory(solstep)
                        history(solstep).q[i][k] -= periodicity.q[k]
                    end
                end
            end
        end
    end
end
