
# Variational Partitioned Runge-Kutta Methods

abstract type VPRKMethod <: VIMethod end

GeometricBase.order(method::VPRKMethod) = RungeKutta.order(tableau(method))

@inline nstages(method::VPRKMethod) = nstages(tableau(method))
@inline eachstage(method::VPRKMethod) = eachstage(tableau(method))

isexplicit(method::VPRKMethod) = RungeKutta.isexplicit(tableau(method))
isimplicit(method::VPRKMethod) = RungeKutta.isimplicit(tableau(method))
issymmetric(method::VPRKMethod) = RungeKutta.issymmetric(tableau(method))
issymplectic(method::VPRKMethod) = RungeKutta.issymplectic(tableau(method))

print_reference(io, method::VPRKMethod) = try ismissing(reference(tableau(method))) || print(io, reference(tableau(method))) catch MethodError String("") end

function Base.show(io::IO, method::VPRKMethod)
    print(io, "\nVariational Partitioned Runge-Kutta Method with Tableau: $(description(tableau(method)))\n")
    print(io, string(tableau(method).q))
    print(io, string(tableau(method).p))
    print_reference(io, method)
end


"""
Variational Partitioned Runge-Kutta Method

```
VPRK(tableau::PartitionedTableau, d=nothing)
VPRK(tableau::Tableau, args...; kwargs...) = VPRK(PartitionedTableau(tableau), args...; kwargs...)
VPRK(tableau1::Tableau, tableau2::Tableau, args...; kwargs...) = VPRK(PartitionedTableau(tableau1, tableau2), args...; kwargs...)
```
"""
struct VPRK{TT,DT} <: VPRKMethod
    tableau::TT
    d::DT

    function VPRK(tableau::TT, d::DT=nothing) where {TT <: PartitionedTableau, DT <: Union{AbstractVector,Nothing}}
        new{TT,DT}(tableau, d)
    end
end

VPRK(tableau::Tableau, args...; kwargs...) = VPRK(PartitionedTableau(tableau), args...; kwargs...)
VPRK(tableau1::Tableau, tableau2::Tableau, args...; kwargs...) = VPRK(PartitionedTableau(tableau1, tableau2), args...; kwargs...)

VPRK(method::RKMethod, args...; kwargs...) = VPRK(PartitionedTableau(tableau(method)))
VPRK(method::PRKMethod, args...; kwargs...) = VPRK(tableau(method))
VPRK(method::VPRKMethod, args...; kwargs...) = VPRK(tableau(method), nullvector(method))

Base.hash(method::VPRK, h::UInt) = hash(method.tableau, hash(method.d, hash(:VPRK, h)))

Base.:(==)(method1::VPRK, method2::VPRK) = (method1.tableau == method2.tableau && method1.d == method2.d)

GeometricBase.tableau(method::VPRK) = method.tableau
nullvector(method::VPRK) = method.d
GeometricBase.order(method::VPRK) = RungeKutta.order(tableau(method))

hasnullvector(method::VPRK{DT,Nothing}) where {DT} = false
hasnullvector(method::VPRK{DT,<:AbstractVector}) where {DT} = true



@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
the symmetric [`TableauSRK3`](@ref) coefficients with 3 stages
for both $a_{ij}$ and $\bar{a}_{ij}$.
"""
struct VPSRK3 <: VPRKMethod end

GeometricBase.tableau(::VPSRK3) = SymplecticPartitionedTableau(TableauSRK3())


@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
Gauss-Legendre coefficients for both $a_{ij}$ and $\bar{a}_{ij}$.
"""
struct VPRKGauss <: VPRKMethod
    s::Int
end

@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
[`TableauLobattoIII`](@ref) for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
struct VPRKLobattoIII <: VPRKMethod
    s::Int
end

@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
[`TableauLobattoIIIA`](@ref) for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
struct VPRKLobattoIIIA <: VPRKMethod
    s::Int
end

@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
[`TableauLobattoIIIB`](@ref) for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
struct VPRKLobattoIIIB <: VPRKMethod
    s::Int
end

@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
[`TableauLobattoIIIC`](@ref) for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
struct VPRKLobattoIIIC <: VPRKMethod
    s::Int
end

@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
[`TableauLobattoIIID`](@ref) for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
struct VPRKLobattoIIID <: VPRKMethod
    s::Int
end

@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
[`TableauLobattoIIIE`](@ref) for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
struct VPRKLobattoIIIE <: VPRKMethod
    s::Int
end

@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
[`TableauLobattoIIIF`](@ref) for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
struct VPRKLobattoIIIF <: VPRKMethod
    s::Int
end

@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
[`TableauLobattoIIIF̄`](@ref) for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
struct VPRKLobattoIIIF̄ <: VPRKMethod
    s::Int
end

@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
[`TableauLobattoIIIG`](@ref) for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
struct VPRKLobattoIIIG <: VPRKMethod
    s::Int
end

@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
[`TableauRadauIIA`](@ref) for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
struct VPRKRadauIIA <: VPRKMethod
    s::Int
end

@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
[`TableauRadauIIB`](@ref) for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
struct VPRKRadauIIB <: VPRKMethod
    s::Int
end

@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
[`TableauLobattoIIIA`](@ref) for the coefficients $a_{ij}$ and [`TableauLobattoIIIB`](@ref) for the coefficients $\bar{a}_{ij}$.
"""
struct VPRKLobattoIIIAIIIB <: VPRKMethod
    s::Int
end

@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
[`TableauLobattoIIIB`](@ref) for the coefficients $a_{ij}$ and [`TableauLobattoIIIA`](@ref) for the coefficients $\bar{a}_{ij}$.
"""
struct VPRKLobattoIIIBIIIA <: VPRKMethod
    s::Int
end

@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
[`TableauLobattoIIIA`](@ref) for the coefficients $a_{ij}$ and its symplectic conjugate for the coefficients $\bar{a}_{ij}$.
"""
struct VPRKLobattoIIIAIIIĀ <: VPRKMethod
    s::Int
end

@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
[`TableauLobattoIIIB`](@ref) for the coefficients $a_{ij}$ and its symplectic conjugate for the coefficients $\bar{a}_{ij}$.
"""
struct VPRKLobattoIIIBIIIB̄ <: VPRKMethod
    s::Int
end

@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
[`TableauLobattoIIIC`](@ref) for the coefficients $a_{ij}$ and its symplectic conjugate for the coefficients $\bar{a}_{ij}$.
"""
struct VPRKLobattoIIICIIIC̄ <: VPRKMethod
    s::Int
end

@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
[`TableauLobattoIIIC̄`](@ref) for the coefficients $a_{ij}$ and its symplectic conjugate for the coefficients $\bar{a}_{ij}$.
"""
struct VPRKLobattoIIIC̄IIIC <: VPRKMethod
    s::Int
end

@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
[`TableauLobattoIIID`](@ref) for the coefficients $a_{ij}$ and its symplectic conjugate for the coefficients $\bar{a}_{ij}$.
"""
struct VPRKLobattoIIIDIIID̄ <: VPRKMethod
    s::Int
end

@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
[`TableauLobattoIIIE`](@ref) for the coefficients $a_{ij}$ and its symplectic conjugate for the coefficients $\bar{a}_{ij}$.
"""
struct VPRKLobattoIIIEIIIĒ <: VPRKMethod
    s::Int
end

@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
[`TableauLobattoIIIF`](@ref) for the coefficients $a_{ij}$ and its symplectic conjugate for the coefficients $\bar{a}_{ij}$.
"""
struct VPRKLobattoIIIFIIIF̄ <: VPRKMethod
    s::Int
end

@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
[`TableauLobattoIIIF̄`](@ref) for the coefficients $a_{ij}$ and its symplectic conjugate for the coefficients $\bar{a}_{ij}$.
"""
struct VPRKLobattoIIIF̄IIIF <: VPRKMethod
    s::Int
end

@doc raw"""
Variational Partitioned Runge-Kutta Method that uses 
[`TableauLobattoIIIG`](@ref) for the coefficients $a_{ij}$ and its symplectic conjugate for the coefficients $\bar{a}_{ij}$.
"""
struct VPRKLobattoIIIGIIIḠ <: VPRKMethod
    s::Int
end


GeometricBase.tableau(method::VPRKGauss) = SymplecticPartitionedTableau(TableauGauss(method.s))
GeometricBase.tableau(method::VPRKLobattoIII) = PartitionedTableau(TableauLobattoIII(method.s))
GeometricBase.tableau(method::VPRKLobattoIIIA) = PartitionedTableau(TableauLobattoIIIA(method.s))
GeometricBase.tableau(method::VPRKLobattoIIIB) = PartitionedTableau(TableauLobattoIIIB(method.s))
GeometricBase.tableau(method::VPRKLobattoIIIC) = PartitionedTableau(TableauLobattoIIIC(method.s))
GeometricBase.tableau(method::VPRKLobattoIIID) = PartitionedTableau(TableauLobattoIIID(method.s))
GeometricBase.tableau(method::VPRKLobattoIIIE) = PartitionedTableau(TableauLobattoIIIE(method.s))
GeometricBase.tableau(method::VPRKLobattoIIIF) = PartitionedTableau(TableauLobattoIIIF(method.s))
GeometricBase.tableau(method::VPRKLobattoIIIF̄) = PartitionedTableau(TableauLobattoIIIF̄(method.s))
GeometricBase.tableau(method::VPRKLobattoIIIG) = PartitionedTableau(TableauLobattoIIIG(method.s))
GeometricBase.tableau(method::VPRKRadauIIA) = PartitionedTableau(TableauRadauIIA(method.s))
GeometricBase.tableau(method::VPRKRadauIIB) = PartitionedTableau(TableauRadauIIB(method.s))

GeometricBase.tableau(method::VPRKLobattoIIIAIIIB) = TableauLobattoIIIAIIIB(method.s)
GeometricBase.tableau(method::VPRKLobattoIIIBIIIA) = TableauLobattoIIIBIIIA(method.s)
GeometricBase.tableau(method::VPRKLobattoIIIAIIIĀ) = TableauLobattoIIIAIIIĀ(method.s)
GeometricBase.tableau(method::VPRKLobattoIIIBIIIB̄) = TableauLobattoIIIBIIIB̄(method.s)
GeometricBase.tableau(method::VPRKLobattoIIICIIIC̄) = TableauLobattoIIICIIIC̄(method.s)
GeometricBase.tableau(method::VPRKLobattoIIIC̄IIIC) = TableauLobattoIIIC̄IIIC(method.s)
GeometricBase.tableau(method::VPRKLobattoIIIDIIID̄) = TableauLobattoIIIDIIID̄(method.s)
GeometricBase.tableau(method::VPRKLobattoIIIEIIIĒ) = TableauLobattoIIIEIIIĒ(method.s)
GeometricBase.tableau(method::VPRKLobattoIIIFIIIF̄) = TableauLobattoIIIFIIIF̄(method.s)
GeometricBase.tableau(method::VPRKLobattoIIIF̄IIIF) = TableauLobattoIIIF̄IIIF(method.s)
GeometricBase.tableau(method::VPRKLobattoIIIGIIIḠ) = TableauLobattoIIIGIIIḠ(method.s)

nullvector(method::VPRKLobattoIIIAIIIB) = get_lobatto_nullvector(method.s)
nullvector(method::VPRKLobattoIIIBIIIA) = get_lobatto_nullvector(method.s)
nullvector(method::VPRKLobattoIIIAIIIĀ) = get_lobatto_nullvector(method.s)
nullvector(method::VPRKLobattoIIIBIIIB̄) = get_lobatto_nullvector(method.s)

GeometricBase.order(::Type{VPRKGauss}) = "2s"
GeometricBase.order(::Type{VPRKLobattoIII}) = "2s-2"
GeometricBase.order(::Type{VPRKLobattoIIIA}) = "2s-2"
GeometricBase.order(::Type{VPRKLobattoIIIB}) = "2s-2"
GeometricBase.order(::Type{VPRKLobattoIIIC}) = "2s-2"
GeometricBase.order(::Type{VPRKLobattoIIID}) = "2s-2"
GeometricBase.order(::Type{VPRKLobattoIIIE}) = "2s-2"
GeometricBase.order(::Type{VPRKLobattoIIIF}) = "2s"
GeometricBase.order(::Type{VPRKLobattoIIIF̄}) = "2s"
GeometricBase.order(::Type{VPRKLobattoIIIG}) = "2s"
GeometricBase.order(::Type{VPRKRadauIIA}) = "2s-1"
GeometricBase.order(::Type{VPRKRadauIIB}) = "2s-1"

GeometricBase.order(::Type{VPRKLobattoIIIAIIIB}) = "2s-2"
GeometricBase.order(::Type{VPRKLobattoIIIBIIIA}) = "2s-2"
GeometricBase.order(::Type{VPRKLobattoIIIAIIIĀ}) = "2s-2"
GeometricBase.order(::Type{VPRKLobattoIIIBIIIB̄}) = "2s-2"
GeometricBase.order(::Type{VPRKLobattoIIICIIIC̄}) = "2s-2"
GeometricBase.order(::Type{VPRKLobattoIIIC̄IIIC}) = "2s-2"
GeometricBase.order(::Type{VPRKLobattoIIIDIIID̄}) = "2s-2"
GeometricBase.order(::Type{VPRKLobattoIIIEIIIĒ}) = "2s-2"
GeometricBase.order(::Type{VPRKLobattoIIIFIIIF̄}) = "2s"
GeometricBase.order(::Type{VPRKLobattoIIIF̄IIIF}) = "2s"
GeometricBase.order(::Type{VPRKLobattoIIIGIIIḠ}) = "2s"

issymplectic(::Type{VPRKGauss}) = true
issymplectic(::Type{VPRKLobattoIII}) = false
issymplectic(::Type{VPRKLobattoIIIA}) = false
issymplectic(::Type{VPRKLobattoIIIB}) = false
issymplectic(::Type{VPRKLobattoIIIC}) = false
issymplectic(::Type{VPRKLobattoIIID}) = true
issymplectic(::Type{VPRKLobattoIIIE}) = true
issymplectic(::Type{VPRKLobattoIIIF}) = false
issymplectic(::Type{VPRKLobattoIIIF̄}) = false
issymplectic(::Type{VPRKLobattoIIIG}) = true
issymplectic(::Type{VPRKRadauIIA}) = false
issymplectic(::Type{VPRKRadauIIB}) = false

issymplectic(::Type{VPRKLobattoIIIAIIIB}) = true
issymplectic(::Type{VPRKLobattoIIIBIIIA}) = true
issymplectic(::Type{VPRKLobattoIIIAIIIĀ}) = true
issymplectic(::Type{VPRKLobattoIIIBIIIB̄}) = true
issymplectic(::Type{VPRKLobattoIIICIIIC̄}) = true
issymplectic(::Type{VPRKLobattoIIIC̄IIIC}) = true
issymplectic(::Type{VPRKLobattoIIIDIIID̄}) = true
issymplectic(::Type{VPRKLobattoIIIEIIIĒ}) = true
issymplectic(::Type{VPRKLobattoIIIFIIIF̄}) = true
issymplectic(::Type{VPRKLobattoIIIF̄IIIF}) = true
issymplectic(::Type{VPRKLobattoIIIGIIIḠ}) = true
