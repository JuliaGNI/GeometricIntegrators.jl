
abstract type AbstractIntegratorRK{dType, tType} <: ODEIntegrator{dType, tType} end
abstract type AbstractIntegratorIRK{dType, tType} <: IODEIntegrator{dType, tType} end
abstract type AbstractIntegratorPRK{dType, tType} <: PODEIntegrator{dType, tType} end

IntegratorRK = Union{AbstractIntegratorRK, AbstractIntegratorIRK, AbstractIntegratorPRK}

@inline GeometricBase.parameters(integrator::IntegratorRK) = integrator.params
@inline GeometricBase.equation(integrator::IntegratorRK, i::Symbol) = parameters(integrator).equs[i]
@inline GeometricBase.equations(integrator::IntegratorRK) = parameters(integrator).equs
@inline GeometricBase.timestep(integrator::IntegratorRK) = parameters(integrator).Δt
@inline tableau(integrator::IntegratorRK)  = parameters(integrator).tab

@inline nstages(integrator::IntegratorRK)  = nstages(tableau(integrator))
@inline eachstage(integrator::IntegratorRK) = 1:nstages(integrator)
