
"Additive Runge Kutta integrator."
immutable IntegratorARK{T} <: Integrator{T}
    equation::Equation
    tableau::TableauARK

end

"Integrate DAE with Additive Runge Kutta integrator."
function integrate!(int::IntegratorARK, s::SolutionDAE)
    # TODO
end
