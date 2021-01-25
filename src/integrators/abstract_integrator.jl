
abstract type Integrator{dType, tType} end

abstract type DeterministicIntegrator{dType, tType} <: Integrator{dType, tType} end

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

parameters(integrator::Integrator) = error("parameters() not implemented for ", typeof(integrator))
equations(integrator::Integrator) = error("equations() not implemented for ", typeof(integrator))
equation(integrator::Integrator) = error("equation() not implemented for ", typeof(integrator))
equation(integrator::Integrator, i::Int) = error("equation() not implemented for ", typeof(integrator))
timestep(integrator::Integrator) = error("timestep() not implemented for ", typeof(integrator))
Base.ndims(integrator::Integrator) = error("ndims() not implemented for ", typeof(integrator))
Common.nconstraints(integrator::Integrator) = error("nconstraints() not implemented for ", typeof(integrator))
nstages(integrator::Integrator) = error("nstages() not implemented for ", typeof(integrator))

eachdim(integrator::Integrator) = 1:ndims(integrator)

get_internal_variables(::Integrator) = NamedTuple()
get_internal_variables(::Nothing) = NamedTuple()


"Create AtomicSolution for ODE."
function Solutions.AtomicSolution(solution::SolutionODE, integrator::Integrator)
    AtomicSolutionODE(get_initial_conditions(solution, 1)..., get_internal_variables(integrator))
end

"Create AtomicSolution for partitioned ODE."
function Solutions.AtomicSolution(solution::SolutionPODE, integrator::Integrator)
    AtomicSolutionPODE(get_initial_conditions(solution, 1)..., get_internal_variables(integrator))
end

"Create AtomicSolution for DAE."
function Solutions.AtomicSolution(solution::SolutionDAE, integrator::Integrator)
    AtomicSolutionDAE(get_initial_conditions(solution, 1)..., get_internal_variables(integrator))
end

"Create AtomicSolution for partitioned DAE."
function Solutions.AtomicSolution(solution::SolutionPDAE, integrator::Integrator)
    AtomicSolutionPDAE(get_initial_conditions(solution, 1)..., get_internal_variables(integrator))
end


abstract type Parameters{DT,TT} end

function_stages!(::Vector{DT}, ::Vector{DT}, ::PT) where {DT, TT, PT <: Parameters{DT,TT}} = error("function_stages!() not implemented for ", PT)
solution_stages!(::Vector{DT}, ::Vector{DT}, ::PT) where {DT, TT, PT <: Parameters{DT,TT}} = error("solution_stages!() not implemented for ", PT)

initialize!(::Integrator, ::AtomicSolution) = nothing

integrate_step!(integrator::Integrator, ::AtomicSolution) = error("integrate_step()! not implemented for ", typeof(integrator))
