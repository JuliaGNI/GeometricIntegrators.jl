"""
`GeometricMethod` is the abstract supertype for all integration methods implemented in GeometricIntegrators.
"""
abstract type GeometricMethod <: AbstractMethod end

abstract type ODEMethod <: GeometricMethod end
abstract type PODEMethod <: GeometricMethod end
abstract type HODEMethod <: GeometricMethod end
abstract type IODEMethod <: GeometricMethod end
abstract type LODEMethod <: GeometricMethod end
abstract type SODEMethod <: GeometricMethod end

abstract type DAEMethod <: GeometricMethod end
abstract type PDAEMethod <: GeometricMethod end
abstract type HDAEMethod <: GeometricMethod end
abstract type IDAEMethod <: GeometricMethod end
abstract type LDAEMethod <: GeometricMethod end

internal_variables(::GeometricMethod, ::GeometricProblem) = NamedTuple()
nullvector(::GeometricMethod) = nothing
GeometricBase.tableau(::GeometricMethod) = missing

default_solver(::GeometricMethod) = NoSolver()
default_iguess(::GeometricMethod) = NoInitialGuess()
default_projection(::GeometricMethod) = NoProjection()

initmethod(method::GeometricMethod) = method
initmethod(method::GeometricMethod, ::AbstractProblem) = initmethod(method)

isodemethod(::Union{GeometricMethod, Type{<:GeometricMethod}}) = false
ispodemethod(::Union{GeometricMethod, Type{<:GeometricMethod}}) = false
ishodemethod(::Union{GeometricMethod, Type{<:GeometricMethod}}) = false
isiodemethod(::Union{GeometricMethod, Type{<:GeometricMethod}}) = false
islodemethod(::Union{GeometricMethod, Type{<:GeometricMethod}}) = false
issodemethod(::Union{GeometricMethod, Type{<:GeometricMethod}}) = false

isdaemethod(::Union{GeometricMethod, Type{<:GeometricMethod}}) = false
ispdaemethod(::Union{GeometricMethod, Type{<:GeometricMethod}}) = false
ishdaemethod(::Union{GeometricMethod, Type{<:GeometricMethod}}) = false
isidaemethod(::Union{GeometricMethod, Type{<:GeometricMethod}}) = false
isldaemethod(::Union{GeometricMethod, Type{<:GeometricMethod}}) = false

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

isexplicit(::GeometricMethod) = missing
isimplicit(::GeometricMethod) = missing
issymmetric(::GeometricMethod) = missing
issymplectic(::GeometricMethod) = missing
isenergypreserving(::GeometricMethod) = missing
isstifflyaccurate(::GeometricMethod) = missing

isexplicit(t::Type{<:GeometricMethod}) = applicable(t) ? isexplicit(t()) : missing
isimplicit(t::Type{<:GeometricMethod}) = applicable(t) ? isimplicit(t()) : missing
issymmetric(t::Type{<:GeometricMethod}) = applicable(t) ? issymmetric(t()) : missing
issymplectic(t::Type{<:GeometricMethod}) = applicable(t) ? issymplectic(t()) : missing
isenergypreserving(t::Type{<:GeometricMethod}) = applicable(t) ? isenergypreserving(t()) : missing
isstifflyaccurate(t::Type{<:GeometricMethod}) = applicable(t) ? isstifflyaccurate(t()) : missing

print_reference(io, method::GeometricMethod) = try ismissing(reference(method)) || print(io, reference(method)) catch MethodError String("") end

# function check_symplecticity end
function symplecticity_conditions end
