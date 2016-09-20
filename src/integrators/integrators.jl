
abstract Integrator{T}

"Integrator: Create integrator for explicit Runge-Kutta tableau."
function Integrator(equation::ODE, tableau::TableauERK, Δt)
    IntegratorERK(equation, tableau, Δt)
end

"Integrator: Create integrator for diagonally implicit Runge-Kutta tableau."
function Integrator(equation::ODE, tableau::TableauDIRK, Δt)
    IntegratorDIRK(equation, tableau, Δt)
end

"Integrator: Create integrator for fully implicit Runge-Kutta tableau."
function Integrator(equation::ODE, tableau::TableauFIRK, Δt)
    IntegratorFIRK(equation, tableau, Δt)
end

"Integrator: Create integrator for partitioned Runge-Kutta tableau."
function Integrator(equation::PODE, tableau::TableauPRK, Δt)
    IntegratorPRK(equation, tableau, Δt)
end

"Integrator: Create integrator for special additive Runge-Kutta tableau."
function Integrator(equation::DAE, tableau::TableauSARK, Δt)
    IntegratorSARK(equation, tableau, Δt)
end

"Integrator: Create integrator for special partitioned additive Runge-Kutta tableau."
function Integrator(equation::PDAE, tableau::TableauSPARK, Δt)
    IntegratorSPARK(equation, tableau, Δt)
end

"Integrator: Print error for integrators not implemented, yet."
function Integrator(equation::Equation, tableau::Tableau, Δt)
    error("No integrator found for tableau ", tableau)
end

"solve: Solve given equation with given tableau for ntime time steps and return solution."
function integrate(equation::Equation, tableau::Tableau, Δt, ntime::Int, nsave::Int=1)
    return integrate(Integrator(equation, tableau, Δt), ntime, nsave)
end

"solve: Apply integrator for ntime time steps and return solution."
function integrate(integrator::Integrator, ntime::Int, nsave::Int=1)
    solution = Solution(integrator.equation, ntime, nsave)
    integrate!(integrator, solution)
    return solution
end


# TODO Add solver status information to all integrators.
