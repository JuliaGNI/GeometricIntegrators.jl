
abstract type IntegratorRK{dType, tType} <: DeterministicIntegrator{dType, tType} end
abstract type IntegratorPRK{dType, tType} <: IntegratorRK{dType, tType} end

@inline equation(integrator::IntegratorRK) = integrator.params.equ
@inline timestep(integrator::IntegratorRK) = integrator.params.Δt
@inline tableau(integrator::IntegratorRK)  = integrator.params.tab

@inline nstages(integrator::IntegratorRK)  = nstages(tableau(integrator))
@inline eachstage(integrator::IntegratorRK) = 1:nstages(integrator)
