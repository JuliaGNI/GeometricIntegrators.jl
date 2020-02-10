
abstract type AbstractIntegratorSPARK{dType, tType} <: IntegratorPRK{dType, tType} end
abstract type AbstractIntegratorHSPARK{dType, tType} <: AbstractIntegratorSPARK{dType, tType} end
abstract type AbstractIntegratorVSPARK{dType, tType} <: AbstractIntegratorSPARK{dType, tType} end
