@doc raw"""
`TableauEPRK`: Tableau of an Explicit Partitioned Runge-Kutta method
```math
\begin{align*}
V_{n,i} &= \hphantom{-} \dfrac{\partial H}{\partial p} (Q_{n,i}, P_{n,i}) , &
Q_{n,i} &= q_{n} + h \sum \limits_{j=1}^{s} a_{ij} \, V_{n,j} , &
q_{n+1} &= q_{n} + h \sum \limits_{i=1}^{s} b_{i} \, V_{n,i} , \\
F_{n,i} &= - \dfrac{\partial H}{\partial q} (Q_{n,i}, P_{n,i}) , &
P_{n,i} &= p_{n} + h  \sum \limits_{i=1}^{s} \bar{a}_{ij} \, F_{n,j} , &
p_{n+1} &= p_{n} + h \sum \limits_{i=1}^{s} \bar{b}_{i} \, F_{n,i} ,
\end{align*}
```
usually satisfying the symplecticity conditions
```math
\begin{align*}
b_{i} \bar{a}_{ij} + b_{j} a_{ji} &= b_{i} b_{j} , &
\bar{b}_i &= b_i .
\end{align*}
```
"""
struct TableauEPRK{T} <: AbstractTableauPRK{T}
    @HeaderTableau

    q::CoefficientsRK{T}
    p::CoefficientsRK{T}

    function TableauEPRK{T}(name, o, q, p) where {T}
        @assert q.s==p.s
        # TODO check that both tableaus are lower triangular and that only one element
        #      a_q[i,i] or a_p[i,i] is non-zero for all i.
        new(name, o, q.s, q, p)
    end
end

function TableauEPRK(name::Symbol, order::Int, q::CoefficientsRK{T}, p::CoefficientsRK{T}) where {T}
    TableauEPRK{T}(name, order, q, p)
end

function TableauEPRK(name::Symbol, order::Int, q::CoefficientsRK{T}) where {T}
    TableauEPRK{T}(name, order, q, q)
end

# TODO function readAbstractTableauEPRKFromFile(dir::AbstractString, name::AbstractString)


"Parameters for right-hand side function of explicit partitioned Runge-Kutta methods."
struct ParametersEPRK{DT, TT, ET <: NamedTuple, D, S} <: Parameters{DT,TT}
    equs::ET
    tab::TableauEPRK{TT}
    Δt::TT

    function ParametersEPRK{DT,D}(equs::ET, tab::TableauEPRK{TT}, Δt::TT) where {D, DT, TT, ET <: NamedTuple}
        new{DT, TT, ET, D, tab.s}(equs, tab, Δt)
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


"Explicit partitioned Runge-Kutta integrator."
struct IntegratorEPRK{DT, TT, PT <: ParametersEPRK{DT,TT}, D, S} <: IntegratorPRK{DT,TT}
    params::PT
    cache::IntegratorCacheEPRK{DT,D,S}

    function IntegratorEPRK{DT,D}(v::Function, f::Function, tableau::TableauEPRK{TT}, Δt::TT) where {DT,TT,D}
        # get number of stages
        S = tableau.s

        # create equation tuple
        equs = NamedTuple{(:v,:f)}((v,f))

        # create params
        params = ParametersEPRK{DT,D}(equs, tableau, Δt)

        # create cache
        cache = IntegratorCacheEPRK{DT,D,S}()

        # create integrators
        new{DT, TT, typeof(params), D, S}(params, cache)
    end

    function IntegratorEPRK(equation::PODE{DT,TT}, tableau::TableauEPRK{TT}, Δt::TT) where {DT,TT}
        IntegratorEPRK{DT, equation.d}(equation.v, equation.f, tableau, Δt)
    end
end


@inline Base.ndims(int::IntegratorEPRK{DT,TT,PT,D,S}) where {DT,TT,PT,D,S} = D


"Compute Q stages of explicit partitioned Runge-Kutta methods."
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
    equations(int)[:f](t, cache.Q[i], cache.P[jmax], cache.F[i])
end

"Compute P stages of explicit partitioned Runge-Kutta methods."
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
    equations(int)[:v](t, cache.Q[jmax], cache.P[i], cache.V[i])
end

"Integrate partitioned ODE with explicit partitioned Runge-Kutta integrator."
function integrate_step!(int::IntegratorEPRK{DT,TT}, sol::AtomicSolutionPODE{DT,TT}) where {DT,TT}
    local tqᵢ::TT
    local tpᵢ::TT

    # store previous solution
    int.cache.Q[0] .= sol.q
    int.cache.P[0] .= sol.p

    # reset cache
    reset!(sol, timestep(int))

    # compute internal stages
    for i in eachstage(int)
        tqᵢ = sol.t̅ + timestep(int) * tableau(int).q.c[i]
        tpᵢ = sol.t̅ + timestep(int) * tableau(int).p.c[i]

        if tableau(int).q.a[i,i] ≠ zero(TT) && tableau(int).p.a[i,i] ≠ zero(TT)
            error("This is an implicit method!")
        elseif tableau(int).q.a[i,i] ≠ zero(TT)
            computeStageP!(int, int.cache, i, i-1, tpᵢ)
            computeStageQ!(int, int.cache, i, i, tqᵢ)
        elseif tableau(int).p.a[i,i] ≠ zero(TT)
            computeStageQ!(int, int.cache, i, i-1, tqᵢ)
            computeStageP!(int, int.cache, i, i, tpᵢ)
        else
            computeStageQ!(int, int.cache, i, i-1, tqᵢ)
            computeStageP!(int, int.cache, i, i-1, tpᵢ)
        end
    end

    # compute final update
    update_solution!(sol.q, sol.q̃, int.cache.V, tableau(int).q.b, tableau(int).q.b̂, timestep(int))
    update_solution!(sol.p, sol.p̃, int.cache.F, tableau(int).p.b, tableau(int).p.b̂, timestep(int))
end
