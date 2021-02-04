
abstract type AbstractIntegratorRK{dType, tType} <: ODEIntegrator{dType, tType} end
abstract type AbstractIntegratorIRK{dType, tType} <: IODEIntegrator{dType, tType} end
abstract type AbstractIntegratorPRK{dType, tType} <: PODEIntegrator{dType, tType} end

IntegratorRK = Union{AbstractIntegratorRK, AbstractIntegratorIRK, AbstractIntegratorPRK}

@inline equation(integrator::IntegratorRK, i::Symbol) = integrator.params.equs[i]
@inline equations(integrator::IntegratorRK) = integrator.params.equs
@inline timestep(integrator::IntegratorRK) = integrator.params.Δt
@inline tableau(integrator::IntegratorRK)  = integrator.params.tab

@inline nstages(integrator::IntegratorRK)  = nstages(tableau(integrator))
@inline eachstage(integrator::IntegratorRK) = 1:nstages(integrator)
