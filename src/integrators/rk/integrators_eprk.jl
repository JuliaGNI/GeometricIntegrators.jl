
"Parameters for right-hand side function of explicit partitioned Runge-Kutta methods."
struct ParametersEPRK{DT, TT, D, S, ET <: NamedTuple} <: Parameters{DT,TT}
    equs::ET
    tab::PartitionedTableau{TT}
    Δt::TT

    function ParametersEPRK{DT,D}(equs::ET, tab::PartitionedTableau{TT}, Δt::TT) where {DT, TT, D, ET <: NamedTuple}
        new{DT, TT, D, tab.s, ET}(equs, tab, Δt)
    end
end


"Explicit Runge-Kutta integrator cache."
struct IntegratorCacheEPRK{DT,D,S} <: PODEIntegratorCache{DT,D}
    Q::OffsetArray{Array{DT,1},1,Array{Array{DT,1},1}}
    P::OffsetArray{Array{DT,1},1,Array{Array{DT,1},1}}

    Y::Vector{Vector{DT}}
    Z::Vector{Vector{DT}}
    V::Vector{Vector{DT}}
    F::Vector{Vector{DT}}

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

function IntegratorCache{ST}(params::ParametersEPRK{DT,TT,D,S}; kwargs...) where {ST,DT,TT,D,S}
    IntegratorCacheEPRK{ST,D,S}(; kwargs...)
end

@inline CacheType(ST, params::ParametersEPRK{DT,TT,D,S}) where {DT,TT,D,S} = IntegratorCacheEPRK{ST,D,S}


@doc raw"""
Explicit partitioned Runge-Kutta integrator solving the system
```math
\begin{aligned}
V_{n,i} &= v (Q_{n,i}, P_{n,i}) , &
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, V_{n,j} , &
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i} \, V_{n,i} , \\
F_{n,i} &= f (Q_{n,i}, P_{n,i}) , &
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
struct IntegratorEPRK{DT, TT, D, S, ET <: NamedTuple} <: AbstractIntegratorPRK{DT,TT}
    params::ParametersEPRK{DT, TT, D, S, ET}
    caches::CacheDict{ParametersEPRK{DT, TT, D, S, ET}}

    function IntegratorEPRK(params::ParametersEPRK{DT,TT,D,S,ET}, caches) where {DT,TT,D,S,ET}
        new{DT, TT, D, S, ET}(params, caches)
    end

    function IntegratorEPRK{DT,D}(equations::ET, tableau::PartitionedTableau{TT}, Δt::TT) where {DT, TT, D, ET <: NamedTuple}
        # get number of stages
        S = tableau.s

        # create params
        params = ParametersEPRK{DT,D}(equations, tableau, Δt)

        # create cache dict
        caches = CacheDict(params)

        # create integrator
        IntegratorEPRK(params, caches)
    end

    function IntegratorEPRK{DT,D}(v::Function, f::Function, tableau::PartitionedTableau{TT}, Δt::TT; kwargs...) where {DT,TT,D}
        IntegratorEPRK{DT,D}(NamedTuple{(:v,:f)}((v,f)), tableau, Δt; kwargs...)
    end

    function IntegratorEPRK{DT,D}(v::Function, f::Function, h::Function, tableau::PartitionedTableau{TT}, Δt::TT; kwargs...) where {DT,TT,D}
        IntegratorEPRK{DT,D}(NamedTuple{(:v,:f,:h)}((v,f,h)), tableau, Δt; kwargs...)
    end

    function IntegratorEPRK(problem::Union{PODEProblem{DT}, HODEProblem{DT}}, tableau::PartitionedTableau{TT}; kwargs...) where {DT,TT}
        IntegratorEPRK{DT, ndims(problem)}(functions(problem), tableau, timestep(problem); kwargs...)
    end
end


@inline Base.ndims(::IntegratorEPRK{DT,TT,D,S}) where {DT,TT,D,S} = D


# Compute Q stages of explicit partitioned Runge-Kutta methods.
function computeStageQ!(int::IntegratorEPRK{DT,TT}, cache::IntegratorCacheEPRK{DT}, i::Int, jmax::Int, t::TT) where {DT,TT}
    fill!(cache.Y[i], zero(DT))
    for j in 1:jmax
        for k in eachdim(int)
            cache.Y[i][k] += tableau(int).q.a[i,j] * cache.V[j][k]
        end
    end
    for k in eachdim(int)
        cache.Q[i][k] = cache.Q[0][k] + timestep(int) * cache.Y[i][k]
    end
    equations(int)[:f](cache.F[i], t, cache.Q[i], cache.P[jmax])
end

# Compute P stages of explicit partitioned Runge-Kutta methods.
function computeStageP!(int::IntegratorEPRK{DT,TT}, cache::IntegratorCacheEPRK{DT}, i::Int, jmax::Int, t::TT) where {DT,TT}
    fill!(cache.Z[i], zero(DT))
    for j in 1:jmax
        for k in eachdim(int)
            cache.Z[i][k] += tableau(int).p.a[i,j] * cache.F[j][k]
        end
    end
    for k in eachdim(int)
        cache.P[i][k] = cache.P[0][k] + timestep(int) * cache.Z[i][k]
    end
    equations(int)[:v](cache.V[i], t, cache.Q[jmax], cache.P[i])
end


function integrate_step!(int::IntegratorEPRK{DT,TT}, sol::SolutionStepPODE{DT,TT},
                         cache::IntegratorCacheEPRK{DT}=int.caches[DT]) where {DT,TT}
    # temporary variables
    local tqᵢ::TT
    local tpᵢ::TT

    # store previous solution
    cache.Q[0] .= sol.q
    cache.P[0] .= sol.p

    # reset atomic solution
    reset!(sol)

    # compute internal stages
    for i in eachstage(int)
        tqᵢ = sol.t̄ + timestep(int) * tableau(int).q.c[i]
        tpᵢ = sol.t̄ + timestep(int) * tableau(int).p.c[i]

        if tableau(int).q.a[i,i] ≠ zero(TT) && tableau(int).p.a[i,i] ≠ zero(TT)
            error("This is a fully implicit method!")
        elseif tableau(int).q.a[i,i] ≠ zero(TT)
            computeStageP!(int, cache, i, i-1, tpᵢ)
            computeStageQ!(int, cache, i, i, tqᵢ)
        elseif tableau(int).p.a[i,i] ≠ zero(TT)
            computeStageQ!(int, cache, i, i-1, tqᵢ)
            computeStageP!(int, cache, i, i, tpᵢ)
        else
            computeStageQ!(int, cache, i, i-1, tqᵢ)
            computeStageP!(int, cache, i, i-1, tpᵢ)
        end
    end

    # compute final update
    update_solution!(sol.q, sol.q̃, cache.V, tableau(int).q.b, tableau(int).q.b̂, timestep(int))
    update_solution!(sol.p, sol.p̃, cache.F, tableau(int).p.b, tableau(int).p.b̂, timestep(int))
end
