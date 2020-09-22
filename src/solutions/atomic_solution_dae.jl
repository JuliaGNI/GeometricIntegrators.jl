"""
Atomic solution for an DAE.

### Fields

* `t`: time of current time step
* `t̅`: time of previous time step
* `q`: current solution of q
* `q̅`: previous solution of q
* `q̃`: compensated summation error of q
* `λ`: current solution of λ
* `λ̅`: previous solution of λ
* `v`: vector field of q
* `v̅`: vector field of q̅
* `u`: projective vector field of q
* `u̅`: projective vector field of q̅
"""
mutable struct AtomicSolutionDAE{DT,TT,IT} <: AtomicSolution{DT,TT}
    t::TT
    t̅::TT

    q::Vector{DT}
    q̅::Vector{DT}
    q̃::Vector{DT}

    λ::Vector{DT}
    λ̅::Vector{DT}

    v::Vector{DT}
    v̅::Vector{DT}

    u::Vector{DT}
    u̅::Vector{DT}

    internal::IT

    function AtomicSolutionDAE{DT, TT, IT}(nd, nm, internal::IT) where {DT <: Number, TT <: Real, IT <: NamedTuple}
        new(zero(TT), zero(TT), zeros(DT, nd), zeros(DT, nd), zeros(DT, nd),
                                zeros(DT, nm), zeros(DT, nm),
                                zeros(DT, nd), zeros(DT, nd),
                                zeros(DT, nd), zeros(DT, nd),
                                internal)
    end
end

AtomicSolutionDAE(DT, TT, nd, nm, internal::IT=NamedTuple()) where {IT} = AtomicSolutionDAE{DT, TT, IT}(nd, nm, internal)

function set_solution!(asol::AtomicSolutionDAE, sol)
    t, q, λ = sol
    asol.t  = t
    asol.q .= q
    asol.λ .= λ
    asol.v .= 0
end

function get_solution(asol::AtomicSolutionDAE)
    (asol.t, asol.q, asol.λ)
end

function CommonFunctions.reset!(asol::AtomicSolutionDAE, Δt)
    asol.t̅  = asol.t
    asol.q̅ .= asol.q
    asol.λ̅ .= asol.λ
    asol.v̅ .= asol.v
    asol.u̅ .= asol.u
    asol.t += Δt
end

function update!(asol::AtomicSolutionDAE{DT}, y::Vector{DT}, λ::Vector{DT}) where {DT}
    for k in eachindex(y)
        update!(asol, y[k], k)
    end
    for k in eachindex(λ)
        asol.λ[k] = λ[k]
    end
end

function update!(asol::AtomicSolutionDAE{DT}, y::DT, k::Int) where {DT}
    asol.q[k], asol.q̃[k] = compensated_summation(y, asol.q[k], asol.q̃[k])
end
