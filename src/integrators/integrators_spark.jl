
"Special Partitioned Additive Runge Kutta integrator."
immutable IntegratorSPARK{T} <: Integrator{T}
    equation::Equation
    tableau::TableauSPARK
    Δt::T

    solver::NonlinearSolver{T}

    x::Array{T,1}
    y::Array{T,1}
    X::Array{T,2}
    Y::Array{T,2}
    F::Array{T,2}
end

"Integrate partitioned DAE with Special Additive Runge Kutta integrator."
function integrate!(int::IntegratorSPARK, s::SolutionPDAE)
    # TODO
end
