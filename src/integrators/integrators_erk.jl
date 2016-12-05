
"Explicit Runge-Kutta integrator."
immutable IntegratorERK{DT,TT,FT} <: Integrator{DT,TT}
    equation::ODE{DT,TT,FT}
    tableau::TableauERK{TT}
    Δt::TT

    x::Array{DT,1}
    y::Array{DT,1}
    X::Array{DT,2}
    Y::Array{DT,2}
    F::Array{DT,2}
    tX::Array{DT,1}
    tF::Array{DT,1}


    function IntegratorERK(equation, tableau, Δt)
        D = equation.d
        S = tableau.s
        new(equation, tableau, Δt,
            zeros(DT,D), zeros(DT,D),
            zeros(DT,D,S), zeros(DT,D,S), zeros(DT,D,S),
            zeros(DT,D), zeros(DT,D))
    end
end

function IntegratorERK{DT,TT,FT}(equation::ODE{DT,TT,FT}, tableau::TableauERK{TT}, Δt::TT)
    IntegratorERK{DT,TT,FT}(equation, tableau, Δt)
end

"Integrate ODE with explicit Runge-Kutta integrator."
function integrate!{DT,TT,FT}(int::IntegratorERK{DT,TT,FT}, sol::SolutionODE{DT})
    local tᵢ::TT

    # loop over initial conditions
    for m in 1:sol.n0
        # copy initial conditions from solution
        simd_copy_xy_first!(int.x, sol, 0, m)

        # loop over time steps
        for n in 1:sol.ntime
            # compute internal stages
            fill!(int.Y, zero(DT))
            for i in 1:int.tableau.s
                for k in 1:sol.nd
                    for j = 1:i-1
                        int.Y[k,i] += int.tableau.a[i,j] * int.F[k,j]
                    end
                    int.X[k,i] = int.x[k] + int.Δt * int.Y[k,i]
                end
                tᵢ = sol.t[n] + int.Δt * int.tableau.c[i]
                simd_copy_xy_first!(int.tX, int.X, i)
                int.equation.f(tᵢ, int.tX, int.tF)
                simd_copy_yx_first!(int.tF, int.F, i)
            end

            # compute final update
            simd_mult!(int.y, int.F, int.tableau.b)
            simd_axpy!(int.Δt, int.y, int.x)

            # copy to solution
            if mod(n, sol.nsave) == 0
                simd_copy_yx_first!(int.x, sol, div(n, sol.nsave), m)
            end
        end
    end
    nothing
end
