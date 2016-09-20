
"Special Additive Runge Kutta integrator."
immutable IntegratorSARK{T} <: Integrator{T}
    equation::Equation
    tableau::TableauSARK

end

"Integrate DAE with Special Additive Runge Kutta integrator."
function integrate!(int::IntegratorSARK, s::SolutionDAE)
    # TODO
end
