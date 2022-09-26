
"""
Abstract atomic or single-step solution.

## Constructors:

```julia
AtomicSolution(equation)
AtomicSolution(solution)
AtomicSolution(solution, integrator)
```

Automatically construct the appropriate atomic solution based on the 
given `equation` or `solution` type. If an `integrator` is provided as,
the `internal` field of the atomic solution is constructed according to
the internal state of the integrator as obtained from the function
`get_internal_variables`.

"""
abstract type AtomicSolution{dType <: Number, tType <: Real, aType <: AbstractArray{dType}} end

"""
Copy the initial conditions of a `GeometricProblem` to the current state of an atomic solution.
"""
function Base.copy!(asol::AtomicSolution, equ::GeometricProblem)
    copy!(asol, initial_conditions(equ))
end

"""
Returns a NamedTuple with the solution of the current time step.
"""
function current end

"""
Returns a NamedTuple with the solution of the previous time step.
"""
function previous end

function GeometricBase.cut_periodic_solution!(asol::AtomicSolution, periodicity)
    if periodicity.q != NullPeriodicity()
        @assert axes(asol.q) == axes(periodicity.q)
        for k in eachindex(asol.q, asol.q̄, periodicity.q)
            if periodicity.q[k] ≠ 0
                while asol.q[k] < 0
                    asol.q[k], asol.q̃[k] = compensated_summation(+periodicity.q[k], asol.q[k], asol.q̃[k])
                    asol.q̄[k] += periodicity.q[k]
                end
                while asol.q[k] ≥ periodicity.q[k]
                    asol.q[k], asol.q̃[k] = compensated_summation(-periodicity.q[k], asol.q[k], asol.q̃[k])
                    asol.q̄[k] -= periodicity.q[k]
                end
            end
        end
    end
end
