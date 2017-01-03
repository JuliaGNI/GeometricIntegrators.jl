
abstract Integrator{dType, tType}

"Create integrator for explicit Runge-Kutta tableau."
function Integrator(equation::ODE, tableau::TableauERK, Δt)
    IntegratorERK(equation, tableau, Δt)
end

"Create integrator for diagonally implicit Runge-Kutta tableau."
function Integrator(equation::ODE, tableau::TableauDIRK, Δt)
    IntegratorDIRK(equation, tableau, Δt)
end

"Create integrator for fully implicit Runge-Kutta tableau."
function Integrator(equation::ODE, tableau::TableauFIRK, Δt)
    IntegratorFIRK(equation, tableau, Δt)
end

"Create integrator for singly implicit Runge-Kutta tableau."
function Integrator(equation::ODE, tableau::TableauSIRK, Δt)
    IntegratorSIRK(equation, tableau, Δt)
end

"Create integrator for explicit partitioned Runge-Kutta tableau."
function Integrator(equation::PODE, tableau::TableauEPRK, Δt)
    IntegratorEPRK(equation, tableau, Δt)
end

"Create integrator for implicit partitioned Runge-Kutta tableau."
function Integrator(equation::PODE, tableau::TableauIPRK, Δt)
    IntegratorIPRK(equation, tableau, Δt)
end

"Create integrator for variational partitioned Runge-Kutta tableau."
function Integrator(equation::IODE, tableau::TableauVPRK, Δt)
    IntegratorVPRK(equation, tableau, Δt)
end

"Create integrator for special additive Runge-Kutta tableau."
function Integrator(equation::DAE, tableau::TableauSARK, Δt)
    IntegratorSARK(equation, tableau, Δt)
end

"Create integrator for special partitioned additive Runge-Kutta tableau."
function Integrator(equation::PDAE, tableau::TableauSPARK, Δt)
    IntegratorSPARK(equation, tableau, Δt)
end

"Create integrator for implicit partitioned additive Runge-Kutta tableau."
function Integrator(equation::IDAE, tableau::TableauIPARK, Δt)
    IntegratorIPARK(equation, tableau, Δt)
end

"Print error for integrators not implemented, yet."
function Integrator(equation::Equation, tableau::AbstractTableau, Δt)
    error("No integrator found for tableau ", tableau)
end

"Apply integrator for ntime time steps and return solution."
function integrate(integrator::Integrator, ntime::Int, nsave::Int=1)
    solution = Solution(integrator.equation, integrator.Δt, ntime, nsave)
    integrate!(integrator, solution)
    return solution
end

"Integrate given equation with given tableau for ntime time steps and return solution."
function integrate(equation::Equation, tableau::AbstractTableau, Δt, ntime::Int, nsave::Int=1)
    return integrate(Integrator(equation, tableau, Δt), ntime, nsave)
end

"Integrate ODE specified by vector field and initial condition with given tableau for ntime time steps and return solution."
function integrate(f::Function, x₀::Vector, tableau::AbstractTableau, Δt, ntime::Int, nsave::Int=1; t₀=0)
    return integrate(ODE(f, t₀, x₀), tableau, Δt, ntime, nsave)
end


# TODO Add solver status information to all integrators (if requested).
