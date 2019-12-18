
abstract type IntegratorCache{DT,D} end

abstract type ODEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type DAEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type IODEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type IDAEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type PODEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type PDAEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end


create_integrator_cache(integrator::Integrator) = error("create_integrator_cache()! not implemented for ", typeof(integrator))

initialize!(::Integrator, ::AtomicSolution) = nothing

integrate_step!(integrator::Integrator, ::AtomicSolution) = error("integrate_step()! not implemented for ", typeof(integrator))

