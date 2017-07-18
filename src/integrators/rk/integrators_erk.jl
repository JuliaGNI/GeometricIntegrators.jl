
"Explicit Runge-Kutta integrator."
struct IntegratorERK{DT,TT,FT} <: Integrator{DT,TT}
    equation::ODE{DT,TT,FT}
    tableau::TableauERK{TT}
    Δt::TT

    q::Array{DT,1}
    V::Array{DT,2}
    tQ::Array{DT,1}
    tV::Array{DT,1}


    function IntegratorERK{DT,TT,FT}(equation, tableau, Δt) where {DT,TT,FT}
        D = equation.d
        S = tableau.q.s
        new(equation, tableau, Δt,
            zeros(DT,D), zeros(DT,D,S),
            zeros(DT,D), zeros(DT,D))
    end
end

function IntegratorERK(equation::ODE{DT,TT,FT}, tableau::TableauERK{TT}, Δt::TT) where {DT,TT,FT}
    IntegratorERK{DT,TT,FT}(equation, tableau, Δt)
end

function initialize!(int::IntegratorERK, sol::SolutionODE, m::Int)
    @assert m ≥ 1
    @assert m ≤ sol.ni

    # copy initial conditions from solution
    get_initial_conditions!(sol, int.q, m)
end

"Integrate ODE with explicit Runge-Kutta integrator."
function integrate_step!(int::IntegratorERK{DT,TT,FT}, sol::SolutionODE{DT,TT,N}, m::Int, n::Int) where {DT,TT,FT,N}
    local tᵢ::TT
    local y::DT

    # compute internal stages
    int.equation.v(sol.t[0] + (n-1)*int.Δt, int.q, int.tV)
    simd_copy_yx_first!(int.tV, int.V, 1)

    for i in 2:int.tableau.q.s
        @inbounds for k in eachindex(int.tQ)
            y = 0
            for j = 1:i-1
                y += int.tableau.q.a[i,j] * int.V[k,j]
            end
            int.tQ[k] = int.q[k] + int.Δt * y
        end
        tᵢ = sol.t[0] + (n-1)*int.Δt + int.Δt * int.tableau.q.c[i]
        int.equation.v(tᵢ, int.tQ, int.tV)
        simd_copy_yx_first!(int.tV, int.V, i)
    end

    # compute final update
    update_solution!(int.q, int.V, int.tableau.q.b, int.Δt)

    # take care of periodic solutions
    cut_periodic_solution!(int.q, int.equation.periodicity)

    # copy to solution
    copy_solution!(sol, int.q, n, m)
end
