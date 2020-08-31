
abstract type Integrator{dType, tType} end

abstract type DeterministicIntegrator{dType, tType} <: Integrator{dType, tType} end
abstract type StochasticIntegrator{dType, tType} <: Integrator{dType, tType} end

abstract type ODEIntegrator{dType, tType} <: DeterministicIntegrator{dType, tType} end
abstract type DAEIntegrator{dType, tType} <: DeterministicIntegrator{dType, tType} end
abstract type PODEIntegrator{dType, tType} <: DeterministicIntegrator{dType, tType} end
abstract type PDAEIntegrator{dType, tType} <: DeterministicIntegrator{dType, tType} end

abstract type IODEIntegrator{dType, tType} <: PODEIntegrator{dType, tType} end
abstract type IDAEIntegrator{dType, tType} <: PDAEIntegrator{dType, tType} end
abstract type HODEIntegrator{dType, tType} <: PODEIntegrator{dType, tType} end
abstract type HDAEIntegrator{dType, tType} <: PDAEIntegrator{dType, tType} end
abstract type VODEIntegrator{dType, tType} <: IODEIntegrator{dType, tType} end
abstract type VDAEIntegrator{dType, tType} <: IDAEIntegrator{dType, tType} end

abstract type SDEIntegrator{dType, tType} <: StochasticIntegrator{dType, tType} end
abstract type PSDEIntegrator{dType, tType} <: StochasticIntegrator{dType, tType} end
abstract type SPSDEIntegrator{dType, tType} <: StochasticIntegrator{dType, tType} end

equation(integrator::Integrator) = error("equation() not implemented for ", typeof(integrator))
timestep(integrator::Integrator) = error("timestep() not implemented for ", typeof(integrator))
Base.ndims(integrator::Integrator) = error("ndims() not implemented for ", typeof(integrator))
CommonFunctions.nconstraints(integrator::Integrator) = error("nconstraints() not implemented for ", typeof(integrator))
noisedims(integrator::Integrator) = error("noisedims() not implemented for ", typeof(integrator))
nstages(integrator::Integrator) = error("nstages() not implemented for ", typeof(integrator))

eachdim(integrator::Integrator) = 1:ndims(integrator)

get_internal_variables(::Integrator) = NamedTuple()


Solutions.AtomicSolution(integrator::ODEIntegrator{DT,TT}) where {DT,TT} =
    AtomicSolutionODE(DT, TT, ndims(integrator), get_internal_variables(integrator))

Solutions.AtomicSolution(integrator::PODEIntegrator{DT,TT}) where {DT,TT} =
    AtomicSolutionPODE(DT, TT, ndims(integrator), get_internal_variables(integrator))

Solutions.AtomicSolution(integrator::DAEIntegrator{DT,TT}) where {DT,TT} =
    AtomicSolutionDAE(DT, TT, ndims(integrator), nconstraints(integrator), get_internal_variables(integrator))

Solutions.AtomicSolution(integrator::PDAEIntegrator{DT,TT}) where {DT,TT} =
    AtomicSolutionPDAE(DT, TT, ndims(integrator), nconstraints(integrator), get_internal_variables(integrator))

Solutions.AtomicSolution(integrator::SDEIntegrator{DT,TT}) where {DT,TT} =
    AtomicSolutionSDE(DT, TT, ndims(integrator), noisedims(integrator), get_internal_variables(integrator))

Solutions.AtomicSolution(integrator::PSDEIntegrator{DT,TT}) where {DT,TT} =
    AtomicSolutionPSDE(DT, TT, ndims(integrator), noisedims(integrator), get_internal_variables(integrator))

Solutions.AtomicSolution(integrator::SPSDEIntegrator{DT,TT}) where {DT,TT} =
    AtomicSolutionPSDE(DT, TT, ndims(integrator), noisedims(integrator), get_internal_variables(integrator))


abstract type Parameters{DT,TT} end

function_stages!(::Vector{DT}, ::Vector{DT}, ::PT) where {DT, TT, PT <: Parameters{DT,TT}} = error("function_stages!() not implemented for ", PT)
solution_stages!(::Vector{DT}, ::Vector{DT}, ::PT) where {DT, TT, PT <: Parameters{DT,TT}} = error("solution_stages!() not implemented for ", PT)

initialize!(::Integrator, ::AtomicSolution) = nothing

integrate_step!(integrator::Integrator, ::AtomicSolution) = error("integrate_step()! not implemented for ", typeof(integrator))
