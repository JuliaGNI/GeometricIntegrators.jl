
abstract type IntegratorCache{DT,D} end

abstract type ODEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type DAEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type SDEIntegratorCache{DT,D,M} <: IntegratorCache{DT,D} end
abstract type IODEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type IDAEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type PODEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type PDAEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type PSDEIntegratorCache{DT,D,M} <: IntegratorCache{DT,D} end


IntegratorCache(params::Parameters) = error("IntegratorCache(params) not implemented for ", typeof(params))
IntegratorCache(integrator::Integrator) = error("IntegratorCache(int)! not implemented for ", typeof(integrator))
