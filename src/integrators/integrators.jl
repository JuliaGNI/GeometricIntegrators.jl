
abstract Integrator

"Integrator: Create integrator appropriate for given tableau."
function Integrator(equation::Equation, tableau::Tableau)
    if typeof(tableau) <: TableauERK
        IntegratorERK(equation, tableau)
    elseif typeof(tableau) <: TableauIRK
        IntegratorIRK(equation, tableau)
    elseif typeof(tableau) <: TableauNLIRK
        IntegratorNLIRK(equation, tableau)
    elseif typeof(tableau) <: TableauPRK
        IntegratorPRK(equation, tableau)
    elseif typeof(tableau) <: TableauSPARK
        IntegratorSPARK(equation, tableau)
    else
        error("No integrator found for tableau ", tableau)
    end
end


# TODO Add integrator specific data structures.

type IntegratorERK <: Integrator
    equation::Equation
    tableau::TableauERK

end


type IntegratorIRK <: Integrator
    equation::Equation
    tableau::TableauIRK

end


type IntegratorNLIRK <: Integrator
    equation::Equation
    tableau::TableauNLIRK

end


type IntegratorPRK <: Integrator
    equation::Equation
    tableau::TableauPRK

end


type IntegratorSPARK <: Integrator
    equation::Equation
    tableau::TableauSPARK

end
