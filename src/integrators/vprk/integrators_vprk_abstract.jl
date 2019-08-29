
abstract type AbstractParametersVPRK{DT,TT,ET,D,S} <: Parameters{DT,TT} end
abstract type AbstractIntegratorVPRK{DT,TT} <: DeterministicIntegrator{DT,TT} end
abstract type AbstractIntegratorCacheVPRK{DT,D,S} <: IODEIntegratorCache{DT,D} end

equation(integrator::AbstractIntegratorVPRK) = integrator.params.equ
timestep(integrator::AbstractIntegratorVPRK) = integrator.params.Δt
tableau(integrator::AbstractIntegratorVPRK) = integrator.params.tab


function update_params!(params::AbstractParametersVPRK, cache::AbstractIntegratorCacheVPRK)
    # set time for nonlinear solver and copy previous solution
    params.t̅  = cache.t
    params.q̅ .= cache.q
    params.p̅ .= cache.p
end
