
abstract type IntegratorRK{dType, tType} <: ODEIntegrator{dType, tType} end
abstract type IntegratorPRK{dType, tType} <: PODEIntegrator{dType, tType} end

AbstractIntegratorRK = Union{IntegratorRK,IntegratorPRK}

@inline equation(integrator::AbstractIntegratorRK, i::Symbol) = integrator.params.equs[i]
@inline equations(integrator::AbstractIntegratorRK) = integrator.params.equs
@inline timestep(integrator::AbstractIntegratorRK) = integrator.params.Δt
@inline tableau(integrator::AbstractIntegratorRK)  = integrator.params.tab

@inline nstages(integrator::AbstractIntegratorRK)  = nstages(tableau(integrator))
@inline eachstage(integrator::AbstractIntegratorRK) = 1:nstages(integrator)
