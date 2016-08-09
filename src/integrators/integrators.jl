
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
    elseif typeof(tableau) <: TableauSARK
        IntegratorSARK(equation, tableau)
    elseif typeof(tableau) <: TableauSPARK
        IntegratorSPARK(equation, tableau)
    else
        error("No integrator found for tableau ", tableau)
    end
end

function solve(equation::Equation, tableau::Tableau, ntime::Int, nsave::Int=1)
    return solve(Integrator(equation, tableau), ntime, nsave)
end

function solve(integrator::Integrator, ntime::Int, nsave::Int=1)
    solution = Solution(integrator, ntime, nsave)
    solve!(integrator, solution)
    return solution
end


# TODO Add integrator specific data structures.

type IntegratorERK <: Integrator
    equation::Equation
    tableau::TableauERK

end

function solve!(int::IntegratorERK, s::SolutionODE)
    # TODO
end

function solve!(int::IntegratorERK, s::SolutionPODE)
    # TODO
end


type IntegratorIRK <: Integrator
    equation::Equation
    tableau::TableauIRK

end

function solve!(int::IntegratorIRK, s::SolutionODE)
    # TODO
end

function solve!(int::IntegratorIRK, s::SolutionPODE)
    # TODO
end


type IntegratorNLIRK <: Integrator
    equation::Equation
    tableau::TableauNLIRK

end

function solve!(int::IntegratorNLIRK, s::SolutionODE)
    # TODO
end

function solve!(int::IntegratorNLIRK, s::SolutionPODE)
    # TODO
end


type IntegratorPRK <: Integrator
    equation::Equation
    tableau::TableauPRK

end

function solve!(int::IntegratorPRK, s::SolutionPODE)
    # TODO
end


type IntegratorSARK <: Integrator
    equation::Equation
    tableau::TableauSARK

end

function solve!(int::IntegratorSARK, s::SolutionDAE)
    # TODO
end


type IntegratorSPARK <: Integrator
    equation::Equation
    tableau::TableauSPARK

end

function solve!(int::IntegratorSPARK, s::SolutionPDAE)
    # TODO
end
