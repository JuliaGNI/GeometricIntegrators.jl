
"Explicit partitioned Runge-Kutta integrator."
struct IntegratorEPRK{DT,TT,VT,FT} <: Integrator{DT,TT}
    equation::PODE{DT,TT,VT,FT}
    tableau::TableauEPRK{TT}
    Δt::TT

    q::Array{DT,1}
    p::Array{DT,1}
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
        S = tableau.s
        new(equation, tableau, Δt,
            zeros(DT,D), zeros(DT,D),
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
function computeStageQ!(int::IntegratorEPRK, i::Int, jmax::Int, t)
    for j in 1:jmax
        for k in 1:int.equation.d
            int.Y[k,i] += int.tableau.q.a[i,j] * int.V[k,j]
        end
    end
    for k in 1:int.equation.d
        int.Q[k,i] = int.q[k] + int.Δt * int.Y[k,i]
    end
    simd_copy_xy_first!(int.tQ, int.Q, i)
    simd_copy_xy_first!(int.tP, int.P, jmax)
    int.equation.f(t, int.tQ, int.tP, int.tF)
    simd_copy_yx_first!(int.tF, int.F, i)
end

"Compute P stages of explicit partitioned Runge-Kutta methods."
function computeStageP!(int::IntegratorEPRK, i::Int, jmax::Int, t)
    for j in 1:jmax
        for k in 1:int.equation.d
            int.Z[k,i] += int.tableau.p.a[i,j] * int.F[k,j]
        end
    end
    for k in 1:int.equation.d
        int.P[k,i] = int.p[k] + int.Δt * int.Z[k,i]
    end
    simd_copy_xy_first!(int.tQ, int.Q, jmax)
    simd_copy_xy_first!(int.tP, int.P, i)
    int.equation.v(t, int.tQ, int.tP, int.tV)
    simd_copy_yx_first!(int.tV, int.V, i)
end

function initialize!(int::IntegratorEPRK, sol::SolutionPODE, m::Int)
    @assert m ≥ 1
    @assert m ≤ sol.ni

    # copy initial conditions from solution
    get_initial_conditions!(sol, int.q, int.p, m)
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
            computeStageP!(int, i, i-1, tpᵢ)
            computeStageQ!(int, i, i, tqᵢ)
        elseif int.tableau.p.a[i,i] ≠ zero(TT)
            computeStageQ!(int, i, i-1, tqᵢ)
            computeStageP!(int, i, i, tpᵢ)
        else
            computeStageQ!(int, i, i-1, tqᵢ)
            computeStageP!(int, i, i-1, tpᵢ)
        end
    end

    # compute final update
    update_solution!(int.q, int.V, int.tableau.q.b, int.Δt)
    update_solution!(int.p, int.F, int.tableau.p.b, int.Δt)

    # copy to solution
    copy_solution!(sol, int.q, int.p, n, m)
end
