"""
Atomic solution for an SDE.

### Fields

* `t`: time of current time step
* `t̅`: time of previous time step
* `q`: current solution of q
* `q̅`: previous solution of q
* `q̃`: compensated summation error of q
* `p`: current solution of p
* `p̅`: previous solution of p
* `p̃`: compensated summation error of p
* `ΔW`: Wiener process driving the stochastic process q
* `ΔZ`: Wiener process driving the stochastic process q
* `K`:  integer parameter defining the truncation of the increments of the Wiener process (for strong solutions)
"""
mutable struct AtomicSolutionPSDE{DT,TT} <: AtomicSolution{DT,TT}
    t::TT
    t̅::TT

    q::Vector{DT}
    q̅::Vector{DT}
    q̃::Vector{DT}

    p::Vector{DT}
    p̅::Vector{DT}
    p̃::Vector{DT}

    ΔW::Vector{DT}
    ΔZ::Vector{DT}
    K::Int

    function AtomicSolutionPSDE{DT, TT}(nd) where {DT <: Number, TT <: Real}
        new(zero(TT), zero(TT), zeros(DT, nd), zeros(DT, nd), zeros(DT, nd),
                                zeros(DT, nd), zeros(DT, nd), zeros(DT, nd),
                                zeros(DT, nd), zeros(DT, nd))
    end
end

AtomicSolutionPSDE(DT, TT, nd) = AtomicSolutionPSDE{DT, TT}(nd)

function set_solution!(asol::AtomicSolutionPSDE, sol)
    t, q, p = sol
    asol.t  = t
    asol.q .= q
    asol.p .= p
end


function set_increments!(asol::AtomicSolutionPSDE, incs)
    ΔW, ΔZ = sol
    asol.ΔW .= ΔW
    asol.ΔZ .= ΔZ
end


function get_solution(asol::AtomicSolutionPSDE)
    (asol.t, asol.q, asol.p)
end

function CommonFunctions.reset!(asol::AtomicSolutionPSDE, Δt)
    asol.t̅  = asol.t
    asol.q̅ .= asol.q
    asol.p̅ .= asol.p
    asol.t += Δt
end

function update!(asol::AtomicSolutionPSDE{DT}, y::Vector{DT}, z::Vector{DT}) where {DT}
    for k in eachindex(y,z)
        update!(asol, y[k], z[k], k)
    end
end

function update!(asol::AtomicSolutionPSDE{DT}, y::DT, z::DT, k::Int) where {DT}
    asol.q[k], asol.q̃[k] = compensated_summation(y, asol.q[k], asol.q̃[k])
    asol.p[k], asol.p̃[k] = compensated_summation(z, asol.p[k], asol.p̃[k])
end
