
"""
`TableauEPRK`: Tableau of an Explicit Partitioned Runge-Kutta method
```math
\\begin{align*}
V_{n,i} &= \\hphantom{-} \\dfrac{\\partial H}{\\partial p} (Q_{n,i}, P_{n,i}) , &
Q_{n,i} &= q_{n} + h \\sum \\limits_{j=1}^{s} a_{ij} \\, V_{n,j} , &
q_{n+1} &= q_{n} + h \\sum \\limits_{i=1}^{s} b_{i} \\, V_{n,i} , \\\\
F_{n,i} &= - \\dfrac{\\partial H}{\\partial q} (Q_{n,i}, P_{n,i}) , &
P_{n,i} &= p_{n} + h  \\sum \\limits_{i=1}^{s} \\bar{a}_{ij} \\, F_{n,j} , &
p_{n+1} &= p_{n} + h \\sum \\limits_{i=1}^{s} \\bar{b}_{i} \\, F_{n,i} ,
\\end{align*}
```
usually satisfying the symplecticity conditions
```math
\\begin{align*}
b_{i} \\bar{a}_{ij} + b_{j} a_{ji} &= b_{i} b_{j} , &
\\bar{b}_i &= b_i .
\\end{align*}
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

    Q::Array{DT,2}
    P::Array{DT,2}
    Y::Array{DT,2}
    Z::Array{DT,2}
    V::Array{DT,2}
    F::Array{DT,2}

    tQ::Array{DT,1}
    tP::Array{DT,1}
    tV::Array{DT,1}
    tF::Array{DT,1}

    function IntegratorEPRK{DT,TT,VT,FT}(equation, tableau, Δt) where {DT,TT,VT,FT}
        D = equation.d
        M = equation.n
        S = tableau.s

        q = Array{Vector{TwicePrecision{DT}}}(undef, M)
        p = Array{Vector{TwicePrecision{DT}}}(undef, M)

        for i in 1:M
            q[i] = zeros(TwicePrecision{DT},D)
            p[i] = zeros(TwicePrecision{DT},D)
        end

        new(equation, tableau, Δt, q, p,
            zeros(DT,D,S), zeros(DT,D,S),
            zeros(DT,D,S), zeros(DT,D,S),
            zeros(DT,D,S), zeros(DT,D,S),
            zeros(DT,D), zeros(DT,D),
            zeros(DT,D), zeros(DT,D))
    end
end

function IntegratorEPRK(equation::PODE{DT,TT,VT,FT}, tableau::TableauEPRK{TT}, Δt::TT) where {DT,TT,VT,FT}
    IntegratorEPRK{DT,TT,VT,FT}(equation, tableau, Δt)
end


"Compute Q stages of explicit partitioned Runge-Kutta methods."
function computeStageQ!(int::IntegratorEPRK, m::Int, i::Int, jmax::Int, t)
    for j in 1:jmax
        for k in 1:int.equation.d
            int.Y[k,i] += int.tableau.q.a[i,j] * int.V[k,j]
        end
    end
    for k in 1:int.equation.d
        int.Q[k,i] = int.q[m][k] + int.Δt * int.Y[k,i]
    end
    simd_copy_xy_first!(int.tQ, int.Q, i)
    jmax == 0 ? int.tP .= int.p[m] : simd_copy_xy_first!(int.tP, int.P, jmax)
    int.equation.f(t, int.tQ, int.tP, int.tF)
    simd_copy_yx_first!(int.tF, int.F, i)
end

"Compute P stages of explicit partitioned Runge-Kutta methods."
function computeStageP!(int::IntegratorEPRK, m::Int, i::Int, jmax::Int, t)
    for j in 1:jmax
        for k in 1:int.equation.d
            int.Z[k,i] += int.tableau.p.a[i,j] * int.F[k,j]
        end
    end
    for k in 1:int.equation.d
        int.P[k,i] = int.p[m][k] + int.Δt * int.Z[k,i]
    end
    jmax == 0 ? int.tQ .= int.q[m] : simd_copy_xy_first!(int.tQ, int.Q, jmax)
    simd_copy_xy_first!(int.tP, int.P, i)
    int.equation.v(t, int.tQ, int.tP, int.tV)
    simd_copy_yx_first!(int.tV, int.V, i)
end

function initialize!(int::IntegratorEPRK, sol::SolutionPODE, m::Int)
    @assert m ≥ 1
    @assert m ≤ sol.ni

    # copy initial conditions from solution
    get_initial_conditions!(sol, int.q[m], int.p[m], m)
end

"Integrate partitioned ODE with explicit partitioned Runge-Kutta integrator."
function integrate_step!(int::IntegratorEPRK{DT,TT,VT,FT}, sol::SolutionPODE{DT,TT,N}, m::Int, n::Int) where {DT,TT,VT,FT,N}
    local j::Int
    local tqᵢ::TT
    local tpᵢ::TT

    # compute internal stages
    fill!(int.Y, zero(DT))
    fill!(int.Z, zero(DT))

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
