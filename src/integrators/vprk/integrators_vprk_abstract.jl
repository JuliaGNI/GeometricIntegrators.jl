
abstract type AbstractParametersVPRK{DT,TT,ET,D,S} <: Parameters{DT,TT} end
abstract type AbstractIntegratorVPRK{DT,TT} <: DeterministicIntegrator{DT,TT} end
abstract type AbstractIntegratorVPRKwProjection{DT,TT} <: AbstractIntegratorVPRK{DT,TT} end

equation(integrator::AbstractIntegratorVPRK) = integrator.params.equ
timestep(integrator::AbstractIntegratorVPRK) = integrator.params.Δt
tableau(integrator::AbstractIntegratorVPRK) = integrator.params.tab
