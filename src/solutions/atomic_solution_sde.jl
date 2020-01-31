"""
Atomic solution for an SDE.

### Fields

* `t`: time of current time step
* `t̅`: time of previous time step
* `q`: current solution of q
* `q̅`: previous solution of q
* `q̃`: compensated summation error of q
* `ΔW`: Wiener process driving the stochastic process q
* `ΔZ`: Wiener process driving the stochastic process q
* `K`:  integer parameter defining the truncation of the increments of the Wiener process (for strong solutions)
"""
mutable struct AtomicSolutionSDE{DT,TT} <: AtomicSolution{DT,TT}
    t::TT
    t̅::TT

    q::Vector{DT}
    q̅::Vector{DT}
    q̃::Vector{DT}

    ΔW::Vector{DT}
    ΔZ::Vector{DT}
    K::Int

    function AtomicSolutionSDE{DT, TT}(nd, nm) where {DT <: Number, TT <: Real}
        new(zero(TT), zero(TT), zeros(DT, nd), zeros(DT, nd), zeros(DT, nd),
                                zeros(DT, nm), zeros(DT, nm))
    end
end

AtomicSolutionSDE(DT, TT, nd, nm) = AtomicSolutionSDE{DT, TT}(nd, nm)

function set_solution!(asol::AtomicSolutionSDE, sol)
    t, q = sol
    asol.t  = t
    asol.q .= q
end

function set_increments!(asol::AtomicSolutionSDE, incs)
    ΔW, ΔZ = incs
    asol.ΔW .= ΔW
    asol.ΔZ .= ΔZ
end

function get_solution(asol::AtomicSolutionSDE)
    (asol.t, asol.q)
end

function get_increments(asol::AtomicSolutionSDE)
    (asol.ΔW, asol.ΔZ)
end

function get_increments!(asol::AtomicSolutionSDE, ΔW, ΔZ)
    ΔW .= asol.ΔW
    ΔZ .= asol.ΔZ
end

function CommonFunctions.reset!(asol::AtomicSolutionSDE, Δt)
    asol.t̅  = asol.t
    asol.q̅ .= asol.q
    asol.t += Δt
end

function update!(asol::AtomicSolutionSDE{DT}, y::Vector{DT}) where {DT}
    for k in eachindex(y)
        update!(asol, y[k], k)
    end
end

function update!(asol::AtomicSolutionSDE{DT}, y::DT, k::Int) where {DT}
    asol.q[k], asol.q̃[k] = compensated_summation(y, asol.q[k], asol.q̃[k])
end
