
"Explicit partitioned Runge-Kutta integrator."
immutable IntegratorEPRK{T} <: Integrator{T}
    equation::PODE{T}
    tableau::TableauEPRK{T}
    Δt::T

    q::Array{T,1}
    p::Array{T,1}
    y::Array{T,1}
    z::Array{T,1}
    Q::Array{T,2}
    P::Array{T,2}
    Y::Array{T,2}
    Z::Array{T,2}
    F::Array{T,2}
    G::Array{T,2}
    tQ::Array{T,1}
    tP::Array{T,1}
    tF::Array{T,1}
    tG::Array{T,1}

    function IntegratorEPRK(equation, tableau, Δt)
        D = equation.d
        S = tableau.s
        new(equation, tableau, Δt,
            zeros(T,D), zeros(T,D),
            zeros(T,D), zeros(T,D),
            zeros(T,D,S), zeros(T,D,S),
            zeros(T,D,S), zeros(T,D,S),
            zeros(T,D,S), zeros(T,D,S),
            zeros(T,D), zeros(T,D),
            zeros(T,D), zeros(T,D))
    end
end

function IntegratorEPRK{T}(equation::PODE{T}, tableau::TableauEPRK{T}, Δt::T)
    IntegratorEPRK{T}(equation, tableau, Δt)
end


"Compute Q stages of explicit partitioned Runge-Kutta methods."
function computeStageQ!(int::IntegratorEPRK, i::Int, jmax::Int, t)
    for j in 1:jmax
        for k in 1:int.equation.d
            int.Y[k,i] += int.tableau.a_q[i,j] * int.F[k,j]
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
            int.Z[k,i] += int.tableau.a_p[i,j] * int.G[k,j]
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
function integrate!{T}(int::IntegratorEPRK{T}, sol::SolutionPODE{T})
    # loop over initial conditions
    for m in 1:sol.n0
        local j::Int
        local tqᵢ::T
        local tpᵢ::T

        # copy initial conditions from solution
        for i in 1:sol.nd
            int.q[i] = sol[1, i, 0, m]
            int.p[i] = sol[2, i, 0, m]
        end

        for n in 1:sol.ntime
            # compute internal stages
            fill!(int.Y, zero(T))
            fill!(int.Z, zero(T))
            for i in 1:int.tableau.s
                tqᵢ = sol.t[n] + int.Δt * int.tableau.c_q[i]
                tpᵢ = sol.t[n] + int.Δt * int.tableau.c_p[i]

                if int.tableau.a_q[i,i] ≠ 0. && int.tableau.a_p[i,i] ≠ 0.
                    error("This is an implicit method!")
                elseif int.tableau.a_q[i,i] ≠ 0.
                    computeStageP!(int, i, i-1, tpᵢ)
                    computeStageQ!(int, i, i, tqᵢ)
                elseif int.tableau.a_p[i,i] ≠ 0.
                    computeStageQ!(int, i, i-1, tqᵢ)
                    computeStageP!(int, i, i, tpᵢ)
                else
                    computeStageQ!(int, i, i-1, tqᵢ)
                    computeStageP!(int, i, i-1, tpᵢ)
                end
            end

            # compute final update
            simd_mult!(int.y, int.F, int.tableau.b_q)
            simd_axpy!(int.Δt, int.y, int.q)

            simd_mult!(int.z, int.G, int.tableau.b_p)
            simd_axpy!(int.Δt, int.z, int.p)

            # copy to solution
            if mod(n, sol.nsave) == 0
                j = div(n, sol.nsave)
                for i in 1:sol.nd
                    sol[1, i, j, m] = int.q[i]
                    sol[2, i, j, m] = int.p[i]
                end
            end
        end
    end
    nothing
end
