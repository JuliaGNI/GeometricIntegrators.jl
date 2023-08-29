
abstract type VIMethod <: LODEMethod end

isiodemethod(::Union{VIMethod, Type{<:VIMethod}}) = true

default_solver(::VIMethod) = Newton()
default_iguess(::VIMethod) = HermiteExtrapolation()
