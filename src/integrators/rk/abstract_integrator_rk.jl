
abstract type IntegratorRK{dType, tType} <: DeterministicIntegrator{dType, tType} end
abstract type IntegratorPRK{dType, tType} <: IntegratorRK{dType, tType} end

@inline equation(integrator::IntegratorRK, i::Symbol) = integrator.params.equs[i]
@inline equations(integrator::IntegratorRK) = integrator.params.equs
@inline timestep(integrator::IntegratorRK) = integrator.params.Δt
@inline tableau(integrator::IntegratorRK)  = integrator.params.tab

@inline nstages(integrator::IntegratorRK)  = nstages(tableau(integrator))
@inline eachstage(integrator::IntegratorRK) = 1:nstages(integrator)
