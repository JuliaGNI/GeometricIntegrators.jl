
abstract type StochasticIntegratorRK{dType, tType, D, M, S} <: StochasticIntegrator{dType, tType} end
abstract type StochasticIntegratorPRK{dType, tType, D, M, S} <: StochasticIntegratorRK{dType, tType, D, M, S} end

@inline parameters(integrator::StochasticIntegratorRK) = integrator.params
@inline equation(integrator::StochasticIntegratorRK, i::Symbol) = integrator.params.equs[i]
@inline equations(integrator::StochasticIntegratorRK) = integrator.params.equs
@inline timestep(integrator::StochasticIntegratorRK) = integrator.params.Δt
@inline tableau(integrator::StochasticIntegratorRK)  = integrator.params.tab

@inline nstages(integrator::StochasticIntegratorRK{DT,TT,D,M,S}) where {DT,TT,D,M,S} = S
@inline noisedims(integrator::StochasticIntegratorRK{DT,TT,D,M,S}) where {DT,TT,D,M,S} = M
@inline Base.ndims(integrator::StochasticIntegratorRK{DT,TT,D,M,S}) where {DT,TT,D,M,S} = D
@inline Base.eltype(integrator::StochasticIntegratorRK{DT}) where {DT} = DT

@inline eachstage(integrator::StochasticIntegratorRK) = 1:nstages(integrator)
@inline eachnoise(integrator::StochasticIntegratorRK) = 1:noisedims(integrator)
@inline eachdim(integrator::StochasticIntegratorRK) = 1:ndims(integrator)
