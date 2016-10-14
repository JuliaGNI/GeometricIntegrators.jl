
"Singly implicit Runge-Kutta integrator."
immutable IntegratorSIRK{T} <: Integrator{T}
    equation::ODE{T}
    tableau::TableauSIRK{T}
    Δt::T

    x::Array{T,1}
    X::Array{T,2}
    Y::Array{T,2}
    F::Array{T,2}

    function IntegratorSIRK(equation, tableau, Δt)
        D = equation.d
        S = tableau.s
        new(equation, tableau, Δt, zeros(T,D), zeros(T,D,S), zeros(T,D,S), zeros(T,D,S))
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
