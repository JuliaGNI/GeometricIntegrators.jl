
abstract type IntegratorRK{dType, tType} <: DeterministicIntegrator{dType, tType} end
abstract type IntegratorPRK{dType, tType} <: IntegratorRK{dType, tType} end

equation(integrator::IntegratorRK) = integrator.params.equ
timestep(integrator::IntegratorRK) = integrator.params.Δt
tableau(integrator::IntegratorRK)  = integrator.params.tab
nstages(integrator::IntegratorRK)  = integrator.params.tab.s
dims(integrator::IntegratorRK)     = integrator.params.equ.d
