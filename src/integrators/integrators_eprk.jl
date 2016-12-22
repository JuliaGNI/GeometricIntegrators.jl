
"Explicit partitioned Runge-Kutta integrator."
immutable IntegratorEPRK{DT,TT,VT,FT} <: Integrator{DT,TT}
    equation::PODE{DT,TT,VT,FT}
    tableau::TableauEPRK{TT}
    Δt::TT

    q::Array{DT,1}
    p::Array{DT,1}
    y::Array{DT,1}
    z::Array{DT,1}
    Q::Array{DT,2}
    P::Array{DT,2}
    Y::Array{DT,2}
    Z::Array{DT,2}
    F::Array{DT,2}
    G::Array{DT,2}
    tQ::Array{DT,1}
    tP::Array{DT,1}
    tF::Array{DT,1}
    tG::Array{DT,1}

    function IntegratorEPRK(equation, tableau, Δt)
        D = equation.d
        S = tableau.s
        new(equation, tableau, Δt,
            zeros(DT,D), zeros(DT,D),
            zeros(DT,D), zeros(DT,D),
            zeros(DT,D,S), zeros(DT,D,S),
            zeros(DT,D,S), zeros(DT,D,S),
            zeros(DT,D,S), zeros(DT,D,S),
            zeros(DT,D), zeros(DT,D),
            zeros(DT,D), zeros(DT,D))
    end
end

function IntegratorEPRK{DT,TT,VT,FT}(equation::PODE{DT,TT,VT,FT}, tableau::TableauEPRK{TT}, Δt::TT)
    IntegratorEPRK{DT,TT,VT,FT}(equation, tableau, Δt)
end


"Compute Q stages of explicit partitioned Runge-Kutta methods."
function computeStageQ!(int::IntegratorEPRK, i::Int, jmax::Int, t)
    for j in 1:jmax
        for k in 1:int.equation.d
            int.Y[k,i] += int.tableau.q.a[i,j] * int.F[k,j]
        end
    end
    for k in 1:int.equation.d
        int.Q[k,i] = int.q[k] + int.Δt * int.Y[k,i]
    end
    simd_copy_xy_first!(int.tQ, int.Q, i)
    int.equation.f(t, int.tQ, int.tG)
    simd_copy_yx_first!(int.tG, int.G, i)
    nothing
end

"Compute P stages of explicit partitioned Runge-Kutta methods."
function computeStageP!(int::IntegratorEPRK, i::Int, jmax::Int, t)
    for j in 1:jmax
        for k in 1:int.equation.d
            int.Z[k,i] += int.tableau.p.a[i,j] * int.G[k,j]
        end
    end
    for k in 1:int.equation.d
        int.P[k,i] = int.p[k] + int.Δt * int.Z[k,i]
    end
    simd_copy_xy_first!(int.tP, int.P, i)
    int.equation.v(t, int.tP, int.tF)
    simd_copy_yx_first!(int.tF, int.F, i)
    nothing
end

"Integrate partitioned ODE with explicit partitioned Runge-Kutta integrator."
function integrate!{DT,TT,VT,FT,N}(int::IntegratorEPRK{DT,TT,VT,FT}, sol::SolutionPODE{DT,TT,N})
    # loop over initial conditions
    for m in 1:sol.ni
        local j::Int
        local tqᵢ::TT
        local tpᵢ::TT

        # copy initial conditions from solution
        get_initial_conditions!(sol, int.q, int.p, m)

        for n in 1:sol.ntime
            # compute internal stages
            fill!(int.Y, zero(DT))
            fill!(int.Z, zero(DT))
            for i in 1:int.tableau.s
                tqᵢ = sol.t[n] + int.Δt * int.tableau.q.c[i]
                tpᵢ = sol.t[n] + int.Δt * int.tableau.p.c[i]

                if int.tableau.q.a[i,i] ≠ 0. && int.tableau.p.a[i,i] ≠ 0.
                    error("This is an implicit method!")
                elseif int.tableau.q.a[i,i] ≠ 0.
                    computeStageP!(int, i, i-1, tpᵢ)
                    computeStageQ!(int, i, i, tqᵢ)
                elseif int.tableau.p.a[i,i] ≠ 0.
                    computeStageQ!(int, i, i-1, tqᵢ)
                    computeStageP!(int, i, i, tpᵢ)
                else
                    computeStageQ!(int, i, i-1, tqᵢ)
                    computeStageP!(int, i, i-1, tpᵢ)
                end
            end

            # compute final update
            simd_mult!(int.y, int.F, int.tableau.q.b)
            simd_axpy!(int.Δt, int.y, int.q)

            simd_mult!(int.z, int.G, int.tableau.p.b)
            simd_axpy!(int.Δt, int.z, int.p)

            # copy to solution
            copy_solution!(sol, int.q, int.p, n, m)
        end
    end
    nothing
end
