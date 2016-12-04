
"Explicit Runge-Kutta integrator."
immutable IntegratorERK{T} <: Integrator{T}
    equation::ODE{T}
    tableau::TableauERK{T}
    Δt::T

    x::Array{T,1}
    y::Array{T,1}
    X::Array{T,2}
    Y::Array{T,2}
    F::Array{T,2}
    tX::Array{T,1}
    tF::Array{T,1}


    function IntegratorERK(equation, tableau, Δt)
        D = equation.d
        S = tableau.s
        new(equation, tableau, Δt, zeros(T,D), zeros(T,D), zeros(T,D,S), zeros(T,D,S), zeros(T,D,S), zeros(T,D), zeros(T,D))
    end
end

function IntegratorERK{T}(equation::Equation{T}, tableau::TableauERK{T}, Δt::T)
    IntegratorERK{T}(equation, tableau, Δt)
end

"Integrate ODE with explicit Runge-Kutta integrator."
function integrate!{T}(int::IntegratorERK{T}, sol::SolutionODE{T})
    local tᵢ::T

    # loop over initial conditions
    for m in 1:sol.n0
        # copy initial conditions from solution
        simd_copy_xy_first!(int.x, sol, 0, m)

        # loop over time steps
        for n in 1:sol.ntime
            # compute internal stages
            fill!(int.Y, zero(T))
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
