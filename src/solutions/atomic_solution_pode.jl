"""
Atomic solution for an PODE.

### Fields

* `t`: time of current time step
* `t̅`: time of previous time step
* `q`: current solution of q
* `q̅`: previous solution of q
* `q̃`: compensated summation error of q
* `p`: current solution of p
* `p̅`: previous solution of p
* `p̃`: compensated summation error of p
* `v`: vector field of q
* `v̅`: vector field of q̅
* `f`: vector field of p
* `f̅`: vector field of p̅
"""
mutable struct AtomicSolutionPODE{DT,TT,IT} <: AtomicSolution{DT,TT}
    t::TT
    t̅::TT

    q::Vector{DT}
    q̅::Vector{DT}
    q̃::Vector{DT}

    p::Vector{DT}
    p̅::Vector{DT}
    p̃::Vector{DT}

    v::Vector{DT}
    v̅::Vector{DT}
    f::Vector{DT}
    f̅::Vector{DT}

    internal::IT

    function AtomicSolutionPODE{DT, TT, IT}(nd, internal::IT) where {DT <: Number, TT <: Real, IT <: NamedTuple}
        new(zero(TT), zero(TT), zeros(DT, nd), zeros(DT, nd), zeros(DT, nd),
                                zeros(DT, nd), zeros(DT, nd), zeros(DT, nd),
                                zeros(DT, nd), zeros(DT, nd), zeros(DT, nd), zeros(DT, nd),
                                internal)
    end
end

AtomicSolutionPODE(DT, TT, nd, internal::IT=NamedTuple()) where {IT} = AtomicSolutionPODE{DT, TT, IT}(nd, internal)

function set_solution!(asol::AtomicSolutionPODE, sol)
    t, q, p = sol
    asol.t  = t
    asol.q .= q
    asol.p .= p
    asol.v .= 0
    asol.f .= 0
end

function get_solution(asol::AtomicSolutionPODE)
    (asol.t, asol.q, asol.p)
end

function CommonFunctions.reset!(asol::AtomicSolutionPODE, Δt)
    asol.t̅  = asol.t
    asol.q̅ .= asol.q
    asol.p̅ .= asol.p
    asol.v̅ .= asol.v
    asol.f̅ .= asol.f
    asol.t += Δt
end

function update!(asol::AtomicSolutionPODE{DT}, y::Vector{DT}, z::Vector{DT}) where {DT}
    for k in eachindex(y,z)
        update!(asol, y[k], z[k], k)
    end
end

function update!(asol::AtomicSolutionPODE{DT}, y::DT, z::DT, k::Int) where {DT}
    asol.q[k], asol.q̃[k] = compensated_summation(y, asol.q[k], asol.q̃[k])
    asol.p[k], asol.p̃[k] = compensated_summation(z, asol.p[k], asol.p̃[k])
end
