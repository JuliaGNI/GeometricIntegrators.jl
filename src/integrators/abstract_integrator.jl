
abstract type AbstractIntegrator{dType, tType} end

abstract type DeterministicIntegrator{dType, tType} <: AbstractIntegrator{dType, tType} end

abstract type ODEIntegrator{dType, tType} <: DeterministicIntegrator{dType, tType} end
abstract type DAEIntegrator{dType, tType} <: DeterministicIntegrator{dType, tType} end
abstract type PODEIntegrator{dType, tType} <: DeterministicIntegrator{dType, tType} end
abstract type PDAEIntegrator{dType, tType} <: DeterministicIntegrator{dType, tType} end

abstract type IODEIntegrator{dType, tType} <: PODEIntegrator{dType, tType} end
abstract type IDAEIntegrator{dType, tType} <: PDAEIntegrator{dType, tType} end
abstract type HODEIntegrator{dType, tType} <: PODEIntegrator{dType, tType} end
abstract type HDAEIntegrator{dType, tType} <: PDAEIntegrator{dType, tType} end
abstract type LODEIntegrator{dType, tType} <: IODEIntegrator{dType, tType} end
abstract type LDAEIntegrator{dType, tType} <: IDAEIntegrator{dType, tType} end

GeometricBase.parameters(integrator::AbstractIntegrator) = error("parameters() not implemented for ", typeof(integrator))
GeometricBase.equations(integrator::AbstractIntegrator) = error("equations() not implemented for ", typeof(integrator))
GeometricBase.equation(integrator::AbstractIntegrator) = error("equation() not implemented for ", typeof(integrator))
GeometricBase.equation(integrator::AbstractIntegrator, i::Int) = error("equation() not implemented for ", typeof(integrator))
GeometricBase.timestep(integrator::AbstractIntegrator) = error("timestep() not implemented for ", typeof(integrator))
Base.ndims(integrator::AbstractIntegrator) = error("ndims() not implemented for ", typeof(integrator))
GeometricBase.nconstraints(integrator::AbstractIntegrator) = error("nconstraints() not implemented for ", typeof(integrator))
# nstages(integrator::AbstractIntegrator) = error("nstages() not implemented for ", typeof(integrator))

eachdim(integrator::AbstractIntegrator) = 1:ndims(integrator)

"""
```julia
get_internal_variables(::Integrator) = NamedTuple()
```
Returns a `NamedTuple` containing all internal variables of an integrator that
shall be stored in an [`SolutionStep`](@ref). If there is no method for a
specific integrator implemented an empty `NamedTuple()` is returned.
"""
get_internal_variables(::AbstractIntegrator) = NamedTuple()
get_internal_variables(::Nothing) = NamedTuple()


# Create SolutionStep with internal variables of integrator.
function Solutions.SolutionStep(solution::GeometricSolution, integrator::AbstractIntegrator)
    SolutionStep(solution, get_internal_variables(integrator))
end


abstract type Parameters{DT,TT} end

function_stages!(::Vector{DT}, ::Vector{DT}, ::PT) where {DT, TT, PT <: Parameters{DT,TT}} = error("function_stages!() not implemented for ", PT)
solution_stages!(::Vector{DT}, ::Vector{DT}, ::PT) where {DT, TT, PT <: Parameters{DT,TT}} = error("solution_stages!() not implemented for ", PT)

initialize!(::AbstractIntegrator, ::SolutionStep) = nothing

cache(::AbstractIntegrator) = missing

"""
Performs one time step with a given integrator.

```julia
integrate_step!(integrator::Integrator, solstep::SolutionStep)
```

The function accepts two arguments: an integrator and an appropriate [`SolutionStep`](@ref),
which contains the state of the system at the beginning and the end of the time step and possibly
additional information like solver output or the solution at internal stages of a Runge-Kutta
method.
"""
integrate_step!(integrator::AbstractIntegrator, ::SolutionStep) = error("integrate_step!() not implemented for ", typeof(integrator))
