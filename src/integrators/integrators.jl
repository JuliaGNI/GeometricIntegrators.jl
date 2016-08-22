
abstract Integrator

"Integrator: Create integrator for explicit Runge-Kutta tableau."
function Integrator(equation::Equation, tableau::TableauERK)
    IntegratorERK(equation, tableau)
end

"Integrator: Create integrator for diagonally implicit Runge-Kutta tableau."
function Integrator(equation::Equation, tableau::TableauIRK)
    IntegratorIRK(equation, tableau)
end

"Integrator: Create integrator for fully implicit Runge-Kutta tableau."
function Integrator(equation::Equation, tableau::TableauNLIRK)
    IntegratorNLIRK(equation, tableau)
end

"Integrator: Create integrator for partitioned Runge-Kutta tableau."
function Integrator(equation::Equation, tableau::TableauPRK)
    IntegratorPRK(equation, tableau)
end

"Integrator: Create integrator for special additive Runge-Kutta tableau."
function Integrator(equation::Equation, tableau::TableauSARK)
    IntegratorSARK(equation, tableau)
end

"Integrator: Create integrator for special partitioned additive Runge-Kutta tableau."
function Integrator(equation::Equation, tableau::TableauSPARK)
    IntegratorSPARK(equation, tableau)
end

"Integrator: Print error for integrators not implemented, yet."
function Integrator(equation::Equation, tableau::Tableau)
    error("No integrator found for tableau ", tableau)
end

"solve: Solve given equation with given tableau for ntime time steps and return solution."
function solve(equation::Equation, tableau::Tableau, ntime::Int, nsave::Int=1)
    return solve(Integrator(equation, tableau), ntime, nsave)
end

"solve: Apply integrator for ntime time steps and return solution."
function solve(integrator::Integrator, ntime::Int, nsave::Int=1)
    solution = Solution(integrator.equation, ntime, nsave)
    solve!(integrator, solution)
end


"IntegratorERK: Explicit Runge-Kutta integrator."
immutable IntegratorERK{T} <: Integrator
    equation::Equation
    tableau::TableauERK

end

function solve!(int::IntegratorERK, s::SolutionODE)
    # TODO
end

function solve!(int::IntegratorERK, s::SolutionPODE)
    # TODO
end


"IntegratorIRK: Diagonally implicit Runge-Kutta integrator."
immutable IntegratorIRK{T} <: Integrator
    equation::Equation
    tableau::TableauIRK

end

function solve!(int::IntegratorIRK, s::SolutionODE)
    # TODO
end

function solve!(int::IntegratorIRK, s::SolutionPODE)
    # TODO
end


"IntegratorIRK: Fully implicit Runge-Kutta integrator."
immutable IntegratorNLIRK{T} <: Integrator
    equation::Equation
    tableau::TableauNLIRK

end

function solve!(int::IntegratorNLIRK, s::SolutionODE)
    # TODO
end

function solve!(int::IntegratorNLIRK, s::SolutionPODE)
    # TODO
end


immutable IntegratorPRK{T} <: Integrator
    equation::Equation
    tableau::TableauPRK

end

function solve!(int::IntegratorPRK, s::SolutionPODE)
    # TODO
end


immutable IntegratorSARK{T} <: Integrator
    equation::Equation
    tableau::TableauSARK

end

function solve!(int::IntegratorSARK, s::SolutionDAE)
    # TODO
end


immutable IntegratorSPARK{T} <: Integrator
    equation::Equation
    tableau::TableauSPARK

end

function solve!(int::IntegratorSPARK, s::SolutionPDAE)
    # TODO
end
