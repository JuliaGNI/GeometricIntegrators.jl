
abstract type Method end

Integrators.integrate(problem::GeometricProblem, method::Method; kwargs...) = integrate(problem, Integrator(problem, method; kwargs...))


abstract type ODEMethod <: Method end
abstract type PODEMethod <: Method end
abstract type HODEMethod <: PODEMethod end
abstract type IODEMethod <: Method end
abstract type LODEMethod <: IODEMethod end
abstract type SODEMethod <: Method end

abstract type DAEMethod <: Method end
abstract type PDAEMethod <: Method end
abstract type HDAEMethod <: PDAEMethod end
abstract type IDAEMethod <: Method end
abstract type LDAEMethod <: IDAEMethod end

isodemethod(::Method) = false
ispodemethod(::Method) = false
ishodemethod(::Method) = false
isiodemethod(::Method) = false
islodemethod(::Method) = false
issodemethod(::Method) = false

isdaemethod(::Method) = false
ispdaemethod(::Method) = false
ishdaemethod(::Method) = false
isidaemethod(::Method) = false
isldaemethod(::Method) = false

isodemethod(::ODEMethod) = true
ispodemethod(::PODEMethod) = true
ishodemethod(::HODEMethod) = true
isiodemethod(::IODEMethod) = true
islodemethod(::LODEMethod) = true
issodemethod(::SODEMethod) = true

isdaemethod(::DAEMethod) = true
ispdaemethod(::PDAEMethod) = true
ishdaemethod(::HDAEMethod) = true
isidaemethod(::IDAEMethod) = true
isldaemethod(::LDAEMethod) = true

isexplicit(::Method) = false
isimplicit(::Method) = false
issymmetric(::Method) = false
issymplectic(::Method) = false
isenergypreserving(::Method) = false
isstifflyaccurate(::Method) = false

order(::Method) = missing
