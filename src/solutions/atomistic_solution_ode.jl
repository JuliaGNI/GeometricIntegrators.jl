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
mutable struct AtomisticSolutionODE{DT,TT} <: AtomisticSolution{DT,TT}
    t::TT
    t̅::TT

    q::Vector{DT}
    q̅::Vector{DT}
    q̃::Vector{DT}

    v::Vector{DT}
    v̅::Vector{DT}

    function AtomisticSolutionODE{DT, TT}(nd) where {DT <: Number, TT <: Real}
        new(zero(TT), zero(TT), zeros(DT, nd), zeros(DT, nd), zeros(DT, nd),
                                zeros(DT, nd), zeros(DT, nd))
    end
end

AtomisticSolutionODE(DT, TT, nd) = AtomisticSolutionODE{DT, TT}(nd)

function CommonFunctions.reset!(asol::AtomisticSolutionODE, Δt)
    asol.t̅  = asol.t
    asol.q̅ .= asol.q
    asol.v̅ .= asol.v
    asol.t += Δt
end

function update!(asol::AtomisticSolutionODE{DT}, v::Vector{DT}) where {DT}
    for k in eachindex(v)
        update!(asol, v[k])
    end
end

function update!(asol::AtomisticSolutionODE{DT}, v::DT, k::Int) where {DT}
    asol.q[k], asol.q̃[k] = compensated_summation(v, asol.q[k], asol.q̃[k])
end
