"""
Atomistic solution for an ODE.

### Fields

* `t`: time of current time step
* `t̅`: time of previous time step
* `q`: current solution of q
* `q̅`: previous solution of q
* `q̃`: compensated summation error of q
* `v`: vector field of q
* `v̅`: vector field of q̅
"""
mutable struct AtomicSolutionODE{DT,TT} <: AtomicSolution{DT,TT}
    t::TT
    t̅::TT

    q::Vector{DT}
    q̅::Vector{DT}
    q̃::Vector{DT}

    v::Vector{DT}
    v̅::Vector{DT}

    function AtomicSolutionODE{DT, TT}(nd) where {DT <: Number, TT <: Real}
        new(zero(TT), zero(TT), zeros(DT, nd), zeros(DT, nd), zeros(DT, nd),
                                zeros(DT, nd), zeros(DT, nd))
    end
end

AtomicSolutionODE(DT, TT, nd) = AtomicSolutionODE{DT, TT}(nd)

function CommonFunctions.set_solution!(asol::AtomicSolutionODE, sol)
    t, q = sol
    asol.t  = t
    asol.q .= q
    asol.v .= 0
end

function CommonFunctions.get_solution(asol::AtomicSolutionODE)
    (asol.t, asol.q)
end

function CommonFunctions.reset!(asol::AtomicSolutionODE, Δt)
    asol.t̅  = asol.t
    asol.q̅ .= asol.q
    asol.v̅ .= asol.v
    asol.t += Δt
end

function update!(asol::AtomicSolutionODE{DT}, v::Vector{DT}) where {DT}
    for k in eachindex(v)
        update!(asol, v[k])
    end
end

function update!(asol::AtomicSolutionODE{DT}, v::DT, k::Int) where {DT}
    asol.q[k], asol.q̃[k] = compensated_summation(v, asol.q[k], asol.q̃[k])
end
