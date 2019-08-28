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

# TODO function readAbstractTableauPRKFromFile(dir::AbstractString, name::AbstractString)


"Explicit partitioned Runge-Kutta integrator."
struct IntegratorEPRK{DT,TT,VT,FT} <: DeterministicIntegrator{DT,TT}
    equation::PODE{DT,TT,VT,FT}
    tableau::TableauEPRK{TT}
    Δt::TT

    function IntegratorEPRK{DT,TT,VT,FT}(equation, tableau, Δt) where {DT,TT,VT,FT}
        new(equation, tableau, Δt)
    end
end

function IntegratorEPRK(equation::PODE{DT,TT,VT,FT}, tableau::TableauEPRK{TT}, Δt::TT) where {DT,TT,VT,FT}
    IntegratorEPRK{DT,TT,VT,FT}(equation, tableau, Δt)
end

equation(int::IntegratorEPRK) = int.equation
timestep(int::IntegratorEPRK) = int.Δt


"Explicit Runge-Kutta integrator cache."
mutable struct IntegratorCacheEPRK{DT,TT,D,S} <: ODEIntegratorCache{DT,D}
    n::Int
    t::TT
    t̅::TT

    q::Vector{TwicePrecision{DT}}
    q̅::Vector{TwicePrecision{DT}}
    p::Vector{TwicePrecision{DT}}
    p̅::Vector{TwicePrecision{DT}}

    Q::OffsetArray{Array{DT,1},1,Array{Array{DT,1},1}}
    P::OffsetArray{Array{DT,1},1,Array{Array{DT,1},1}}

    Y::Vector{Vector{DT}}
    Z::Vector{Vector{DT}}
    V::Vector{Vector{DT}}
    F::Vector{Vector{DT}}

    function IntegratorCacheEPRK{DT,TT,D,S}() where {DT,TT,D,S}
        q = zeros(TwicePrecision{DT}, D)
        q̅ = zeros(TwicePrecision{DT}, D)
        p = zeros(TwicePrecision{DT}, D)
        p̅ = zeros(TwicePrecision{DT}, D)

        Q = create_internal_stage_vector_with_zero(DT, D, S)
        P = create_internal_stage_vector_with_zero(DT, D, S)

        Y = create_internal_stage_vector(DT, D, S)
        Z = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        F = create_internal_stage_vector(DT, D, S)

        new(0, zero(TT), zero(TT), q, q̅, p, p̅, Q, P, Y, Z, V, F)
    end
end

function create_integrator_cache(int::IntegratorEPRK{DT,TT}) where {DT,TT}
    IntegratorCacheEPRK{DT, TT, ndims(equation(int)), int.tableau.s}()
end

function CommonFunctions.reset!(cache::IntegratorCacheEPRK{DT,TT}, Δt::TT) where {DT,TT}
    cache.t̅  = cache.t
    cache.q̅ .= cache.q
    cache.p̅ .= cache.p
    cache.t += Δt
    cache.n += 1
end

function CommonFunctions.get_solution(cache::IntegratorCacheEPRK)
    (cache.t, cache.q, cache.p)
end

function CommonFunctions.set_solution!(cache::IntegratorCacheEPRK, sol, n=0)
    t, q, p = sol
    cache.n  = n
    cache.t  = t
    cache.q .= q
    cache.p .= p
end


"Compute Q stages of explicit partitioned Runge-Kutta methods."
function computeStageQ!(int::IntegratorEPRK{DT,TT}, cache::IntegratorCacheEPRK{DT,TT}, i::Int, jmax::Int, t) where {DT,TT}
    fill!(cache.Y[i], zero(DT))
    for j in 1:jmax
        for k in 1:int.equation.d
            cache.Y[i][k] += int.tableau.q.a[i,j] * cache.V[j][k]
        end
    end
    for k in 1:int.equation.d
        cache.Q[i][k] = cache.q[k] + int.Δt * cache.Y[i][k]
    end
    int.equation.f(t, cache.Q[i], cache.P[jmax], cache.F[i])
end

"Compute P stages of explicit partitioned Runge-Kutta methods."
function computeStageP!(int::IntegratorEPRK{DT,TT}, cache::IntegratorCacheEPRK{DT,TT}, i::Int, jmax::Int, t) where {DT,TT}
    fill!(cache.Z[i], zero(DT))
    for j in 1:jmax
        for k in 1:int.equation.d
            cache.Z[i][k] += int.tableau.p.a[i,j] * cache.F[j][k]
        end
    end
    for k in 1:int.equation.d
        cache.P[i][k] = cache.p[k] + int.Δt * cache.Z[i][k]
    end
    int.equation.v(t, cache.Q[jmax], cache.P[i], cache.V[i])
end

"Integrate partitioned ODE with explicit partitioned Runge-Kutta integrator."
function integrate_step!(int::IntegratorEPRK{DT,TT}, cache::IntegratorCacheEPRK{DT,TT}) where {DT,TT}
    local tqᵢ::TT
    local tpᵢ::TT

    # store previous solution
    cache.Q[0] .= cache.q
    cache.P[0] .= cache.p

    # compute internal stages
    for i in 1:int.tableau.s
        tqᵢ = cache.t̅ + int.Δt * int.tableau.q.c[i]
        tpᵢ = cache.t̅ + int.Δt * int.tableau.p.c[i]

        if int.tableau.q.a[i,i] ≠ zero(TT) && int.tableau.p.a[i,i] ≠ zero(TT)
            error("This is an implicit method!")
        elseif int.tableau.q.a[i,i] ≠ zero(TT)
            computeStageP!(int, cache, i, i-1, tpᵢ)
            computeStageQ!(int, cache, i, i, tqᵢ)
        elseif int.tableau.p.a[i,i] ≠ zero(TT)
            computeStageQ!(int, cache, i, i-1, tqᵢ)
            computeStageP!(int, cache, i, i, tpᵢ)
        else
            computeStageQ!(int, cache, i, i-1, tqᵢ)
            computeStageP!(int, cache, i, i-1, tpᵢ)
        end
    end

    # compute final update
    update_solution!(cache.q, cache.V, int.tableau.q.b, int.Δt)
    update_solution!(cache.p, cache.F, int.tableau.p.b, int.Δt)

    # take care of periodic solutions
    cut_periodic_solution!(cache.q, int.equation.periodicity)
end
