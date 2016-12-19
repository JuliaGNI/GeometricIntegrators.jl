
"Explicit Runge-Kutta integrator."
immutable IntegratorERK{DT,TT,FT} <: Integrator{DT,TT}
    equation::ODE{DT,TT,FT}
    tableau::TableauERK{TT}
    Δt::TT

    x::Array{DT,1}
    F::Array{DT,2}
    tX::Array{DT,1}
    tF::Array{DT,1}


    function IntegratorERK(equation, tableau, Δt)
        D = equation.d
        S = tableau.s
        new(equation, tableau, Δt,
            zeros(DT,D), zeros(DT,S,D),
            zeros(DT,D), zeros(DT,D))
    end
end

function IntegratorERK{DT,TT,FT}(equation::ODE{DT,TT,FT}, tableau::TableauERK{TT}, Δt::TT)
    IntegratorERK{DT,TT,FT}(equation, tableau, Δt)
end

"Integrate ODE with explicit Runge-Kutta integrator."
function integrate!{DT,TT,FT,N}(int::IntegratorERK{DT,TT,FT}, sol::SolutionODE{DT,TT,N})
    local tᵢ::TT
    local y::DT

    # loop over initial conditions
    for m in 1:sol.ni
        # copy initial conditions from solution
        get_data!(sol.q, int.x, 0, m)

        # loop over time steps
        for n in 1:sol.ntime
            # compute internal stages
            int.equation.v(sol.t[n], int.x, int.tF)
            simd_copy_yx_second!(int.tF, int.F, 1)

            for i in 2:int.tableau.s
                @inbounds for k in eachindex(int.tX)
                    y = 0
                    for j = 1:i-1
                        y += int.tableau.a[i,j] * int.F[j,k]
                    end
                    int.tX[k] = int.x[k] + int.Δt * y
                end
                tᵢ = sol.t[n] + int.Δt * int.tableau.c[i]
                int.equation.v(tᵢ, int.tX, int.tF)
                simd_copy_yx_second!(int.tF, int.F, i)
            end

            # compute final update
            simd_abXpy!(int.Δt, int.tableau.b, int.F, int.x)

            # copy to solution
            set_data!(sol.q, int.x, n, m)
        end
    end
    nothing
end
