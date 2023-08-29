
abstract type HPIMethod <: LODEMethod end

isiodemethod(::Union{HPIMethod, Type{<:HPIMethod}}) = true

default_solver(::HPIMethod) = Newton()
default_iguess(::HPIMethod) = HermiteExtrapolation()
