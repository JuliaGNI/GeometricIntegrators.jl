
"Singly implicit Runge-Kutta integrator."
immutable IntegratorSIRK{DT,TT,FT} <: Integrator{DT,TT}
    equation::ODE{DT,TT,FT}
    tableau::TableauSIRK{TT}
    Δt::TT

    x::Array{DT,1}
    X::Array{DT,2}
    Y::Array{DT,2}
    F::Array{DT,2}

    function IntegratorSIRK(equation, tableau, Δt)
        D = equation.d
        S = tableau.s
        new(equation, tableau, Δt, zeros(DT,D), zeros(DT,D,S), zeros(DT,D,S), zeros(DT,D,S))
    end
end

function IntegratorSIRK(equation::Equation, tableau::TableauSIRK, Δt)
    T = eltype(equation.q₀)
    IntegratorSIRK{T}(equation, tableau, Δt)
end

"Integrate ODE with singly implicit Runge-Kutta integrator."
function integrate!(int::IntegratorSIRK, s::SolutionODE)
    # TODO
end

"Integrate partitioned ODE with singly implicit Runge-Kutta integrator."
function integrate!(int::IntegratorSIRK, s::SolutionPODE)
    # TODO
end
