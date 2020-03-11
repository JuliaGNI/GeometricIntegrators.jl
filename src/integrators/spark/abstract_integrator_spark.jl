
abstract type AbstractIntegratorSPARK{dType, tType, D, S, R} <: IntegratorPRK{dType, tType} end
abstract type AbstractIntegratorHSPARK{dType, tType, D, S, R} <: AbstractIntegratorSPARK{dType, tType, D, S, R} end
abstract type AbstractIntegratorVSPARK{dType, tType, D, S, R} <: AbstractIntegratorSPARK{dType, tType, D, S, R} end

@inline parameters(int::AbstractIntegratorSPARK) = int.params
@inline nstages(int::AbstractIntegratorSPARK{DT,TT,D,S,R}) where {DT,TT,D,S,R} = S
@inline pstages(int::AbstractIntegratorSPARK{DT,TT,D,S,R}) where {DT,TT,D,S,R} = R
@inline Base.ndims(int::AbstractIntegratorSPARK{DT,TT,D,S,R}) where {DT,TT,D,S,R} = D
