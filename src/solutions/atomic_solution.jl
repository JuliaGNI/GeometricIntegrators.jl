
"""
Abstract atomistic or single-step solution.

## Constructors:

```julia
AtomicSolution(equation)
AtomicSolution(solution)
AtomicSolution(solution, integrator)
```

Automatically construct the appropriate atomistic solution based on the 
given `equation` or `solution` type. If an `integrator` is provided as,
the `internal` field of the atomic solution is constructed according to
the internal state of the integrator as obtained from the function
`get_internal_variables`.

"""
abstract type AtomicSolution{dType <: Number, tType <: Real, aType <: AbstractArray{dType}} end


"""
```julia
copy_solution!(solution, atomic_solution, n, m)
```

Copy solution for time step `n` and initial condition `m` from atomic solution
to solution object.

"""
function copy_solution!(sol::AbstractSolution, asol::AtomicSolution, n, m)
    copy_solution!(sol, get_solution(asol)..., n, m)
end

function GeometricBase.cut_periodic_solution!(asol::AtomicSolution{DT,TT,AT}, periodicity::AT) where {DT,TT,AT}
    @assert axes(asol.q) == axes(periodicity)
    for k in eachindex(asol.q, asol.q̄, periodicity)
        if periodicity[k] ≠ 0
            while asol.q[k] < 0
                asol.q[k], asol.q̃[k] = compensated_summation(+periodicity[k], asol.q[k], asol.q̃[k])
                asol.q̄[k] += periodicity[k]
            end
            while asol.q[k] ≥ periodicity[k]
                asol.q[k], asol.q̃[k] = compensated_summation(-periodicity[k], asol.q[k], asol.q̃[k])
                asol.q̄[k] -= periodicity[k]
            end
        end
    end
end
