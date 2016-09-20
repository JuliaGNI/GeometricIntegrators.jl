
immutable IntegratorSARK{T} <: Integrator{T}
    equation::Equation
    tableau::TableauSARK

end

"solve!: Solve DAE with Special Additive Runge Kutta integrator."
function integrate!(int::IntegratorSARK, s::SolutionDAE)
    # TODO
end
