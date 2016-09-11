
immutable IntegratorSARK{T} <: Integrator{T}
    equation::Equation
    tableau::TableauSARK

end

function solve!(int::IntegratorSARK, s::SolutionDAE)
    # TODO
end
