
abstract type AbstractIntegratorRK{dType, tType} <: ODEIntegrator{dType, tType} end
abstract type AbstractIntegratorIRK{dType, tType} <: IODEIntegrator{dType, tType} end
abstract type AbstractIntegratorPRK{dType, tType} <: PODEIntegrator{dType, tType} end

IntegratorRK = Union{AbstractIntegratorRK, AbstractIntegratorIRK, AbstractIntegratorPRK}

# @inline nstages(integrator::IntegratorRK)  = nstages(tableau(integrator))
# @inline eachstage(integrator::IntegratorRK) = 1:nstages(integrator)


GeometricBase.tableau(int::Integrator{<:GeometricProblem, <:RKMethod}) = tableau(method(int))

initmethod(method::RKMethod) = RK(method)
initmethod(method::PRKMethod) = PRK(method)

Integrator(problem::Union{PODEProblem,HODEProblem}, method::RKMethod, args...; kwargs...) = Integrator(problem, PRK(method), args...; kwargs...)


function solver(solver::Type{<:NewtonMethod}, problem::GeometricProblem, method::GeometricMethod, caches::CacheDict; F::Callable = function_stages!)
    # create solution vector for nonlinear solver
    x = zeros(DT, ndims(problem) * nstages(method))

    # create wrapper function f!(x,b) that calls `F(x, b, params)`
    # with the appropriate `params`
    f! = (b,x) -> F(b, x, solstep, problem, method, caches)

    # create nonlinear solver with solver type obtained from config dictionary
    solver(x, zero(x), f!)
end

function update!(solstep::SolutionStepODE, V, tableau::Tableau, Δt)
    update_solution!(solstep.q, solstep.q̄[1], solstep.q̃, V, tableau.b, tableau.b̂, Δt)
end

function update!(solstep::SolutionStepPODE, V, F, tableau::PartitionedTableau, Δt)
    update_solution!(solstep.q, solstep.q̄[1], solstep.q̃, V, tableau.q.b, tableau.q.b̂, Δt)
    update_solution!(solstep.p, solstep.p̄[1], solstep.p̃, F, tableau.p.b, tableau.p.b̂, Δt)
end
