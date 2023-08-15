
"Explicit Runge-Kutta integrator cache."
struct IntegratorCacheEPRK{DT,D,S} <: PODEIntegratorCache{DT,D}
    Q::OffsetArray{Array{DT,1},1,Array{Array{DT,1},1}}
    P::OffsetArray{Array{DT,1},1,Array{Array{DT,1},1}}

    V::Vector{Vector{DT}}
    F::Vector{Vector{DT}}
    Y::Vector{Vector{DT}}
    Z::Vector{Vector{DT}}

    function IntegratorCacheEPRK{DT,D,S}() where {DT,D,S}
        Q = create_internal_stage_vector_with_zero(DT, D, S)
        P = create_internal_stage_vector_with_zero(DT, D, S)

        Y = create_internal_stage_vector(DT, D, S)
        Z = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        F = create_internal_stage_vector(DT, D, S)

        new(Q, P, Y, Z, V, F)
    end
end

function Cache{ST}(problem::EquationProblem, method::EPRK; kwargs...) where {ST}
    S = nstages(tableau(method))
    D = ndims(problem)
    IntegratorCacheEPRK{ST,D,S}(; kwargs...)
end

@inline CacheType(ST, problem::EquationProblem, method::EPRK) = IntegratorCacheEPRK{ST, ndims(problem), nstages(tableau(method))}


@doc raw"""
Explicit partitioned Runge-Kutta integrator solving the system
```math
\begin{aligned}
V_{n,i} &= v (t_i, Q_{n,i}, P_{n,i}) , &
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, V_{n,j} , &
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i} \, V_{n,i} , \\
F_{n,i} &= f (t_i, Q_{n,i}, P_{n,i}) , &
P_{n,i} &= p_{n} + h  \sum \limits_{i=1}^{s} \bar{a}_{ij} \, F_{n,j} , &
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} \bar{b}_{i} \, F_{n,i} .
\end{aligned}
```
Usually we are interested in Hamiltonian systems, where
```math
\begin{aligned}
V_{n,i} &= \dfrac{\partial H}{\partial p} (Q_{n,i}, P_{n,i}) , &
F_{n,i} &= - \dfrac{\partial H}{\partial q} (Q_{n,i}, P_{n,i}) , 
\end{aligned}
```
and tableaus satisfying the symplecticity conditions
```math
\begin{aligned}
b_{i} \bar{a}_{ij} + \bar{b}_{j} a_{ji} &= b_{i} \bar{b}_{j} , &
\bar{b}_i &= b_i .
\end{aligned}
```
"""
const IntegratorEPRK{DT,TT} = Integrator{<:PODEProblem{DT,TT}, <:EPRK}

initmethod(method::EPRKMethod) = EPRK(method)

function Base.show(io::IO, int::IntegratorEPRK)
    print(io, "\nExplicit Partitioned Runge-Kutta Integrator with:\n")
    print(io, "   Timestep: $(timestep(int))\n")
    print(io, "   Tableau:  $(description(tableau(int)))\n")
    print(io, "   $(string(tableau(int).q))")
    print(io, "   $(string(tableau(int).p))")
    # print(io, reference(int.params.tab))
end


# Compute Q stages of explicit partitioned Runge-Kutta methods.
function compute_stage_q!(solstep::SolutionStepPODE{DT,TT}, problem::PODEProblem{DT,TT}, method::EPRK, caches::CacheDict, i::Int, jmax::Int, t::TT) where {DT,TT}
    # obtain cache
    local Q = caches[DT].Q
    local P = caches[DT].P
    local V = caches[DT].V
    local F = caches[DT].F
    local Y = caches[DT].Y

    fill!(Y[i], zero(DT))
    for k in eachindex(Y[i])
        for j in 1:jmax
            Y[i][k] += timestep(problem) * tableau(method).q.a[i,j] * V[j][k]
        end
    end
    for k in eachindex(Y[i])
        Q[i][k] = solstep.q̄[k] + Y[i][k]
    end
    functions(problem).f(F[i], t, Q[i], P[jmax])
end

# Compute P stages of explicit partitioned Runge-Kutta methods.
function compute_stage_p!(solstep::SolutionStepPODE{DT,TT}, problem::PODEProblem{DT,TT}, method::EPRK, caches::CacheDict, i::Int, jmax::Int, t::TT) where {DT,TT}
    # obtain cache
    local Q = caches[DT].Q
    local P = caches[DT].P
    local V = caches[DT].V
    local F = caches[DT].F
    local Z = caches[DT].Z

    fill!(Z[i], zero(DT))
    for k in eachindex(Z[i])
        for j in 1:jmax
            Z[i][k] += timestep(problem) * tableau(method).p.a[i,j] * F[j][k]
        end
    end
    for k in eachindex(Z[i])
        P[i][k] = solstep.p̄[k] + Z[i][k]
    end
    functions(problem).v(V[i], t, Q[jmax], P[i])
end


function integrate_step!(solstep::SolutionStepPODE{DT,TT}, problem::PODEProblem{DT,TT}, method::EPRK, caches::CacheDict, ::NoSolver) where {DT,TT}
    # temporary variables
    local tqᵢ::TT
    local tpᵢ::TT

    # store previous solution
    caches[DT].Q[0] .= solstep.q
    caches[DT].P[0] .= solstep.p

    # compute internal stages
    for i in eachstage(method)
        tqᵢ = solstep.t̄ + timestep(problem) * tableau(method).q.c[i]
        tpᵢ = solstep.t̄ + timestep(problem) * tableau(method).p.c[i]

        if tableau(method).q.a[i,i] ≠ zero(TT) && tableau(method).p.a[i,i] ≠ zero(TT)
            error("This is an implicit method!")
        elseif tableau(method).q.a[i,i] ≠ zero(TT)
            compute_stage_p!(solstep, problem, method, caches, i, i-1, tpᵢ)
            compute_stage_q!(solstep, problem, method, caches, i, i, tqᵢ)
        elseif tableau(method).p.a[i,i] ≠ zero(TT)
            compute_stage_q!(solstep, problem, method, caches, i, i-1, tqᵢ)
            compute_stage_p!(solstep, problem, method, caches, i, i, tpᵢ)
        else
            compute_stage_q!(solstep, problem, method, caches, i, i-1, tqᵢ)
            compute_stage_p!(solstep, problem, method, caches, i, i-1, tpᵢ)
        end
    end

    # compute final update
    update!(solstep, caches[DT].V, caches[DT].F, tableau(method), timestep(problem))
end
