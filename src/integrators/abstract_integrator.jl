
abstract type Integrator{dType, tType} end

abstract type DeterministicIntegrator{dType, tType} <: Integrator{dType, tType} end
abstract type StochasticIntegrator{dType, tType} <: Integrator{dType, tType} end

abstract type ODEIntegrator{dType, tType} <: DeterministicIntegrator{dType, tType} end
abstract type DAEIntegrator{dType, tType} <: DeterministicIntegrator{dType, tType} end
abstract type IODEIntegrator{dType, tType} <: DeterministicIntegrator{dType, tType} end
abstract type IDAEIntegrator{dType, tType} <: DeterministicIntegrator{dType, tType} end
abstract type PODEIntegrator{dType, tType} <: DeterministicIntegrator{dType, tType} end
abstract type PDAEIntegrator{dType, tType} <: DeterministicIntegrator{dType, tType} end

abstract type HODEIntegrator{dType, tType} <: PODEIntegrator{dType, tType} end
abstract type HDAEIntegrator{dType, tType} <: PDAEIntegrator{dType, tType} end
abstract type VODEIntegrator{dType, tType} <: IODEIntegrator{dType, tType} end
abstract type VDAEIntegrator{dType, tType} <: IDAEIntegrator{dType, tType} end

abstract type SDEIntegrator{dType, tType} <: StochasticIntegrator{dType, tType} end
abstract type PSDEIntegrator{dType, tType} <: StochasticIntegrator{dType, tType} end
abstract type SPSDEIntegrator{dType, tType} <: StochasticIntegrator{dType, tType} end

equation(integrator::Integrator) = error("equation() not implemented for ", typeof(integrator))
timestep(integrator::Integrator) = error("timestep() not implemented for ", typeof(integrator))

Base.ndims(integrator::Integrator) = ndims(equation(integrator))

eachdim(integrator::Integrator) = 1:ndims(integrator)


abstract type Parameters{DT,TT} end

function_stages!(x::Vector{DT}, b::Vector{DT}, params::PT) where {DT, TT, PT <: Parameters{DT,TT}} = error("function_stages!() not implemented for ", PT)
