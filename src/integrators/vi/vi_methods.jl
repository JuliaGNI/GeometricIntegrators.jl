
abstract type VIMethod <: LODEMethod end

default_solver(::VIMethod) = Newton()
default_iguess(::VIMethod) = HermiteExtrapolation()
