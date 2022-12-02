
abstract type Method end

Integrators.integrate(problem::GeometricProblem, method::Method; kwargs...) = integrate(problem, Integrator(problem, method; kwargs...))


abstract type ODEMethod <: Method end
abstract type PODEMethod <: Method end
abstract type HODEMethod <: Method end
abstract type IODEMethod <: Method end
abstract type LODEMethod <: Method end
abstract type SODEMethod <: Method end

abstract type DAEMethod <: Method end
abstract type PDAEMethod <: Method end
abstract type HDAEMethod <: Method end
abstract type IDAEMethod <: Method end
abstract type LDAEMethod <: Method end

isodemethod(::Union{Method, Type{<:Method}}) = false
ispodemethod(::Union{Method, Type{<:Method}}) = false
ishodemethod(::Union{Method, Type{<:Method}}) = false
isiodemethod(::Union{Method, Type{<:Method}}) = false
islodemethod(::Union{Method, Type{<:Method}}) = false
issodemethod(::Union{Method, Type{<:Method}}) = false

isdaemethod(::Union{Method, Type{<:Method}}) = false
ispdaemethod(::Union{Method, Type{<:Method}}) = false
ishdaemethod(::Union{Method, Type{<:Method}}) = false
isidaemethod(::Union{Method, Type{<:Method}}) = false
isldaemethod(::Union{Method, Type{<:Method}}) = false

isodemethod(::Union{ODEMethod, Type{<:ODEMethod}}) = true
ispodemethod(::Union{PODEMethod, Type{<:PODEMethod}}) = true
ishodemethod(::Union{HODEMethod, Type{<:HODEMethod}}) = true
isiodemethod(::Union{IODEMethod, Type{<:IODEMethod}}) = true
islodemethod(::Union{LODEMethod, Type{<:LODEMethod}}) = true
issodemethod(::Union{SODEMethod, Type{<:SODEMethod}}) = true

isdaemethod(::Union{DAEMethod, Type{<:DAEMethod}}) = true
ispdaemethod(::Union{PDAEMethod, Type{<:PDAEMethod}}) = true
ishdaemethod(::Union{HDAEMethod, Type{<:HDAEMethod}}) = true
isidaemethod(::Union{IDAEMethod, Type{<:IDAEMethod}}) = true
isldaemethod(::Union{LDAEMethod, Type{<:LDAEMethod}}) = true

isexplicit(::Method) = missing
isimplicit(::Method) = missing
issymmetric(::Method) = missing
issymplectic(::Method) = missing
isenergypreserving(::Method) = missing
isstifflyaccurate(::Method) = missing

isexplicit(t::Type{<:Method}) = applicable(t) ? isexplicit(t()) : missing
isimplicit(t::Type{<:Method}) = applicable(t) ? isimplicit(t()) : missing
issymmetric(t::Type{<:Method}) = applicable(t) ? issymmetric(t()) : missing
issymplectic(t::Type{<:Method}) = applicable(t) ? issymplectic(t()) : missing
isenergypreserving(t::Type{<:Method}) = applicable(t) ? isenergypreserving(t()) : missing
isstifflyaccurate(t::Type{<:Method}) = applicable(t) ? isstifflyaccurate(t()) : missing

order(::Method) = missing
order(t::Type{<:Method}) = applicable(t) ? order(t()) : missing

RungeKutta.description(::Union{Method, Type{<:Method}}) = missing
