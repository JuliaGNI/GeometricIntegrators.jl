
abstract type AbstractIntegratorVPRK{DT,TT,D,S} <: AbstractIntegratorIRK{DT,TT} end
abstract type AbstractIntegratorVPRKwProjection{DT,TT,D,S} <: AbstractIntegratorVPRK{DT,TT,D,S} end

@inline parameters(integrator::AbstractIntegratorVPRK) = integrator.params

@inline nstages(integrator::AbstractIntegratorVPRK{DT,TT,D,S}) where {DT,TT,D,S} = S
@inline Base.ndims(integrator::AbstractIntegratorVPRK{DT,TT,D,S}) where {DT,TT,D,S} = D

@inline eachstage(integrator::AbstractIntegratorVPRK) = 1:nstages(integrator)
@inline eachdim(integrator::AbstractIntegratorVPRK) = 1:ndims(integrator)
