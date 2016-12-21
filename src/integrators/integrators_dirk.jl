
"Diagonally implicit Runge-Kutta integrator."
immutable IntegratorDIRK{T} <: Integrator{T}
    equation::ODE{T}
    tableau::TableauDIRK{T}
    Δt::T

    x::Array{T,1}
    X::Array{T,2}
    Y::Array{T,2}
    F::Array{T,2}

    function IntegratorDIRK(equation, tableau, Δt)
        D = equation.d
        S = tableau.q.s
        new(equation, tableau, Δt, zeros(T,D), zeros(T,D,S), zeros(T,D,S), zeros(T,D,S))
    end
end

function IntegratorDIRK(equation::Equation, tableau::TableauDIRK, Δt)
    T = eltype(equation.q₀)
    IntegratorDIRK{T}(equation, tableau, Δt)
end

"Integrate ODE with diagonally implicit Runge-Kutta integrator."
function integrate!(int::IntegratorDIRK, s::SolutionODE)
    # TODO
end

"Integrate partitioned ODE with diagonally implicit Runge-Kutta integrator."
function integrate!(int::IntegratorDIRK, s::SolutionPODE)
    # TODO
end
