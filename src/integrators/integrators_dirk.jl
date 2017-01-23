
"Diagonally implicit Runge-Kutta integrator."
immutable IntegratorDIRK{DT,TT,FT} <: Integrator{DT,TT}
    equation::ODE{DT,TT,FT}
    tableau::TableauDIRK{TT}
    Δt::TT

    x::Array{DT,1}
    X::Array{DT,2}
    Y::Array{DT,2}
    F::Array{DT,2}

    function IntegratorDIRK(equation, tableau, Δt)
        D = equation.d
        S = tableau.q.s
        new(equation, tableau, Δt, zeros(DT,D), zeros(DT,D,S), zeros(DT,D,S), zeros(DT,D,S))
    end
end

function IntegratorDIRK(equation::Equation, tableau::TableauDIRK, Δt)
    DT = eltype(equation.q₀)
    TT = typeof(Δt)
    FT = typeof(equation.v)
    IntegratorDIRK{DT,TT,FT}(equation, tableau, Δt)
end

"Integrate ODE with diagonally implicit Runge-Kutta integrator."
function integrate!(int::IntegratorDIRK, s::SolutionODE)
    # TODO
end

"Integrate partitioned ODE with diagonally implicit Runge-Kutta integrator."
function integrate!(int::IntegratorDIRK, s::SolutionPODE)
    # TODO
end
