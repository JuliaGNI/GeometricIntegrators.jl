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
struct IntegratorEPRK{DT,TT,VT,FT} <: Integrator{DT,TT}
    equation::PODE{DT,TT,VT,FT}
    tableau::TableauEPRK{TT}
    Δt::TT

    q::Vector{Vector{TwicePrecision{DT}}}
    p::Vector{Vector{TwicePrecision{DT}}}

    Q::OffsetArray{Array{DT,1},1,Array{Array{DT,1},1}}
    P::OffsetArray{Array{DT,1},1,Array{Array{DT,1},1}}

    Y::Vector{Vector{DT}}
    Z::Vector{Vector{DT}}
    V::Vector{Vector{DT}}
    F::Vector{Vector{DT}}

    function IntegratorEPRK{DT,TT,VT,FT}(equation, tableau, Δt) where {DT,TT,VT,FT}
        D = equation.d
        M = equation.n
        S = tableau.s

        q = create_solution_vector(DT, D, M)
        p = create_solution_vector(DT, D, M)

        Q = create_internal_stage_vector_with_zero(DT, D, S)
        P = create_internal_stage_vector_with_zero(DT, D, S)

        Y = create_internal_stage_vector(DT, D, S)
        Z = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        F = create_internal_stage_vector(DT, D, S)

        new(equation, tableau, Δt, q, p, Q, P, Y, Z, V, F)
    end
end

function IntegratorEPRK(equation::PODE{DT,TT,VT,FT}, tableau::TableauEPRK{TT}, Δt::TT) where {DT,TT,VT,FT}
    IntegratorEPRK{DT,TT,VT,FT}(equation, tableau, Δt)
end


"Compute Q stages of explicit partitioned Runge-Kutta methods."
function computeStageQ!(int::IntegratorEPRK{DT,TT}, m::Int, i::Int, jmax::Int, t) where {DT,TT}
    fill!(int.Y[i], zero(DT))
    for j in 1:jmax
        for k in 1:int.equation.d
            int.Y[i][k] += int.tableau.q.a[i,j] * int.V[j][k]
        end
    end
    for k in 1:int.equation.d
        int.Q[i][k] = int.q[m][k] + int.Δt * int.Y[i][k]
    end
    int.equation.f(t, int.Q[i], int.P[jmax], int.F[i])
end

"Compute P stages of explicit partitioned Runge-Kutta methods."
function computeStageP!(int::IntegratorEPRK{DT,TT}, m::Int, i::Int, jmax::Int, t) where {DT,TT}
    fill!(int.Z[i], zero(DT))
    for j in 1:jmax
        for k in 1:int.equation.d
            int.Z[i][k] += int.tableau.p.a[i,j] * int.F[j][k]
        end
    end
    for k in 1:int.equation.d
        int.P[i][k] = int.p[m][k] + int.Δt * int.Z[i][k]
    end
    int.equation.v(t, int.Q[jmax], int.P[i], int.V[i])
end

function initialize!(int::IntegratorEPRK, sol::SolutionPODE, m::Int)
    @assert m ≥ 1
    @assert m ≤ sol.ni

    # copy initial conditions from solution
    get_initial_conditions!(sol, int.q[m], int.p[m], m)
end

"Integrate partitioned ODE with explicit partitioned Runge-Kutta integrator."
function integrate_step!(int::IntegratorEPRK{DT,TT}, sol::SolutionPODE{DT,TT}, m::Int, n::Int) where {DT,TT}
    local j::Int
    local tqᵢ::TT
    local tpᵢ::TT

    # store previous solution
    int.Q[0] .= int.q[m]
    int.P[0] .= int.p[m]

    # compute internal stages
    for i in 1:int.tableau.s
        tqᵢ = sol.t[n] + int.Δt * int.tableau.q.c[i]
        tpᵢ = sol.t[n] + int.Δt * int.tableau.p.c[i]

        if int.tableau.q.a[i,i] ≠ zero(TT) && int.tableau.p.a[i,i] ≠ zero(TT)
            error("This is an implicit method!")
        elseif int.tableau.q.a[i,i] ≠ zero(TT)
            computeStageP!(int, m, i, i-1, tpᵢ)
            computeStageQ!(int, m, i, i, tqᵢ)
        elseif int.tableau.p.a[i,i] ≠ zero(TT)
            computeStageQ!(int, m, i, i-1, tqᵢ)
            computeStageP!(int, m, i, i, tpᵢ)
        else
            computeStageQ!(int, m, i, i-1, tqᵢ)
            computeStageP!(int, m, i, i-1, tpᵢ)
        end
    end

    # compute final update
    update_solution!(int.q[m], int.V, int.tableau.q.b, int.Δt)
    update_solution!(int.p[m], int.F, int.tableau.p.b, int.Δt)

    # take care of periodic solutions
    cut_periodic_solution!(int.q[m], int.equation.periodicity)

    # copy to solution
    copy_solution!(sol, int.q[m], int.p[m], n, m)
end
