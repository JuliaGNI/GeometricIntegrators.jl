
"""
Abstract atomistic or single-step solution.
"""
abstract type AtomicSolution{dType <: Number, tType <: Real, aType <: AbstractArray{dType}} end


"Copy solution from atomistic solution to solution object."
function copy_solution!(sol::Solution, asol::AtomicSolution, n, m)
    copy_solution!(sol, get_solution(asol)..., n, m)
end

function Common.cut_periodic_solution!(asol::AtomicSolution{DT,TT,AT}, periodicity::AT) where {DT,TT,AT}
    @assert axes(asol.q) == axes(periodicity)
    for k in eachindex(asol.q, asol.q̅, periodicity)
        if periodicity[k] ≠ 0
            while asol.q[k] < 0
                asol.q[k], asol.q̃[k] = compensated_summation(+periodicity[k], asol.q[k], asol.q̃[k])
                asol.q̅[k] += periodicity[k]
            end
            while asol.q[k] ≥ periodicity[k]
                asol.q[k], asol.q̃[k] = compensated_summation(-periodicity[k], asol.q[k], asol.q̃[k])
                asol.q̅[k] -= periodicity[k]
            end
        end
    end
end
