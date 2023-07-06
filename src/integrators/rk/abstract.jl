
abstract type AbstractIntegratorRK{dType, tType} <: ODEIntegrator{dType, tType} end
abstract type AbstractIntegratorIRK{dType, tType} <: IODEIntegrator{dType, tType} end
abstract type AbstractIntegratorPRK{dType, tType} <: PODEIntegrator{dType, tType} end

const IntegratorRK = Union{AbstractIntegratorRK, AbstractIntegratorIRK, AbstractIntegratorPRK}

# @inline nstages(integrator::IntegratorRK)  = nstages(tableau(integrator))
# @inline eachstage(integrator::IntegratorRK) = 1:nstages(integrator)

const StageVector{T} = Vector{<:AbstractVector{T}}

GeometricBase.tableau(int::Integrator{<:GeometricProblem, <:RKMethod}) = tableau(method(int))

initmethod(method::RKMethod) = RK(method)
initmethod(method::PRKMethod) = PRK(method)

Integrator(problem::Union{PODEProblem,HODEProblem}, method::RKMethod, args...; kwargs...) = Integrator(problem, PRK(method), args...; kwargs...)
