
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


function update!(sol::SolutionStepODE, V, tableau::Tableau, Δt)
    update_solution!(sol.q, sol.q̃, V, tableau.b, tableau.b̂, Δt)
    sol.t += Δt
end

function update!(sol::SolutionStepPODE, V, F, tableau::PartitionedTableau, Δt)
    update_solution!(sol.q, sol.q̃, V, tableau.q.b, tableau.q.b̂, Δt)
    update_solution!(sol.p, sol.p̃, F, tableau.p.b, tableau.p.b̂, Δt)
    sol.t += Δt
end
