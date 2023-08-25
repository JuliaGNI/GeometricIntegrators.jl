
abstract type RKMethod <: ODEMethod end
abstract type PRKMethod <: PODEMethod end

abstract type ERKMethod <: RKMethod end
abstract type IRKMethod <: RKMethod end
abstract type DIRKMethod <: IRKMethod end

abstract type EPRKMethod <: PRKMethod end
abstract type IPRKMethod <: PRKMethod end

const RungeKuttaMethod = Union{RKMethod,PRKMethod}

GeometricBase.tableau(method::RKMethod) = error("No tableau for Runge-Kutta method $(typeof(method)) provided")
GeometricBase.tableau(method::PRKMethod) = error("No tableau for partitioned Runge-Kutta method $(typeof(method)) provided")

GeometricBase.order(method::RungeKuttaMethod) = order(tableau(method))

nstages(method::RungeKuttaMethod) = RungeKutta.nstages(tableau(method))
eachstage(method::RungeKuttaMethod) = RungeKutta.eachstage(tableau(method))
coefficients(method::RungeKuttaMethod) = RungeKutta.coefficients(tableau(method))
weights(method::RungeKuttaMethod) = RungeKutta.weights(tableau(method))
nodes(method::RungeKuttaMethod) = RungeKutta.nodes(tableau(method))

print_reference(io, method::RungeKuttaMethod) = try ismissing(reference(tableau(method))) || print(io, reference(tableau(method))) catch MethodError String("") end

ispodemethod(::Union{RKMethod, Type{<:RKMethod}}) = true
ishodemethod(::Union{RKMethod, Type{<:RKMethod}}) = true
isiodemethod(::Union{RKMethod, Type{<:RKMethod}}) = true
islodemethod(::Union{RKMethod, Type{<:RKMethod}}) = true
ishodemethod(::Union{PRKMethod, Type{<:PRKMethod}}) = true
isiodemethod(::Union{PRKMethod, Type{<:PRKMethod}}) = true
islodemethod(::Union{PRKMethod, Type{<:PRKMethod}}) = true

isexplicit(method::RungeKuttaMethod) = RungeKutta.isexplicit(tableau(method))
isimplicit(method::RungeKuttaMethod) = RungeKutta.isimplicit(tableau(method))
issymmetric(method::RungeKuttaMethod) = RungeKutta.issymmetric(tableau(method))
issymplectic(method::RungeKuttaMethod) = RungeKutta.issymplectic(tableau(method))
# isenergypreserving(method::RungeKuttaMethod) = RungeKutta.order(tableau(method))
# isstifflyaccurate(method::RungeKuttaMethod) = RungeKutta.order(tableau(method))

isexplicit(method::ERKMethod) = true
isexplicit(method::IRKMethod) = false
isimplicit(method::ERKMethod) = false
isimplicit(method::IRKMethod) = true

isexplicit(method::EPRKMethod) = true
isexplicit(method::IPRKMethod) = false
isimplicit(method::EPRKMethod) = false
isimplicit(method::IPRKMethod) = true


function Base.show(io::IO, method::RKMethod)
    print(io, "\nRunge-Kutta Method with Tableau: $(description(tableau(method)))\n")
    print(io, string(tableau(method)))
    print_reference(io, method)
end

function Base.show(io::IO, method::PRKMethod)
    print(io, "\nPartitioned Runge-Kutta Method with Tableau: $(description(tableau(method)))\n")
    print(io, string(tableau(method).q))
    print(io, string(tableau(method).p))
    print_reference(io, method)
end


"""
Explicit Runge-Kutta Method

```
ERK(tableau)
```
"""
struct ERK{TT <: Tableau} <: ERKMethod
    tableau::TT
end

"""
Implicit Runge-Kutta Method

```
IRK(tableau)
```
"""
struct IRK{TT <: Tableau, ImplicitUpdate} <: IRKMethod
    tableau::TT

    function IRK(tableau::TT; implicit_update::Bool = false) where {TT <: Tableau}
        new{TT, implicit_update}(tableau)
    end
end

implicit_update(::IRKMethod) = false
implicit_update(::IRK{TT,IU}) where {TT,IU} = IU

"""
Diagonally Implicit Runge-Kutta Method

```
DIRK(tableau)
```
"""
struct DIRK{TT <: Tableau} <: DIRKMethod
    tableau::TT
end


"""
Runge-Kutta Method

```
RK(tableau)
```

Returns an explicit, implicit or diagonally implicit Runge-Kutta method depending on the tableau.
"""
function RK(tableau::Tableau)
    if RungeKutta.isexplicit(tableau)
        # Create method for explicit Runge-Kutta tableau
        return ERK(tableau)
    elseif RungeKutta.isdiagonallyimplicit(tableau)
        # Create method for diagonally implicit Runge-Kutta tableau
        return DIRK(tableau)
    elseif RungeKutta.isfullyimplicit(tableau)
        # Create method for implicit Runge-Kutta tableau
        return IRK(tableau)
    end
end

RK(method::RKMethod, args...; kwargs...) = RK(tableau(method), args...; kwargs...)
ERK(method::RKMethod, args...; kwargs...) = ERK(tableau(method), args...; kwargs...)
IRK(method::RKMethod, args...; kwargs...) = IRK(tableau(method), args...; kwargs...)
DIRK(method::RKMethod, args...; kwargs...) = DIRK(tableau(method), args...; kwargs...)

@inline GeometricBase.tableau(method::Union{ERK,IRK,DIRK}) = method.tableau


"""
Explicit Partitioned Runge-Kutta Method

```
EPRK(tableau)
```
"""
struct EPRK{TT <: PartitionedTableau} <: EPRKMethod
    tableau::TT
end

"""
Implicit Partitioned Runge-Kutta Method

```
IPRK(tableau)
```
"""
struct IPRK{TT <: PartitionedTableau} <: IPRKMethod
    tableau::TT
end

"""
Partitioned Runge-Kutta Method

```
PRK(tableau)
```

Returns an explicit or implicit partitioned Runge-Kutta method depending on the tableau.
"""
function PRK(tableau::PartitionedTableau)
    if RungeKutta.isexplicit(tableau)
        # Create method for explicit partitioned Runge-Kutta tableau
        return EPRK(tableau)
    elseif RungeKutta.isimplicit(tableau)
        # Create method for implicit partitioned Runge-Kutta tableau
        return IPRK(tableau)
    end
end

PRK(tableau::Tableau, args...; kwargs...) = PRK(PartitionedTableau(tableau), args...; kwargs...)
PRK(method::PRKMethod, args...; kwargs...) = PRK(tableau(method), args...; kwargs...)
PRK(method::RKMethod, args...; kwargs...) = PRK(tableau(method), args...; kwargs...)
EPRK(method::PRKMethod, args...; kwargs...) = EPRK(tableau(method), args...; kwargs...)
IPRK(method::PRKMethod, args...; kwargs...) = IPRK(tableau(method), args...; kwargs...)

@inline GeometricBase.tableau(method::Union{EPRK,IPRK}) = method.tableau



# Explicit Runge-Kutta Methods

"""
Explicit Runge-Kutta method with [`TableauForwardEuler`](@ref).

$(reference(Val(:ExplicitEuler)))
"""
struct ForwardEuler <: ERKMethod end

"""
Explicit Runge-Kutta method with [`TableauExplicitEuler`](@ref).

$(reference(Val(:ExplicitEuler)))
"""
struct ExplicitEulerRK <: ERKMethod end

"""
Explicit Runge-Kutta method with [`TableauExplicitMidpoint`](@ref).

$(reference(Val(:ExplicitMidpoint)))
"""
struct ExplicitMidpoint <: ERKMethod end

"""
Explicit Runge-Kutta method with [`TableauHeun2`](@ref).

$(reference(Val(:Heun2)))
"""
struct Heun2 <: ERKMethod end

"""
Explicit Runge-Kutta method with [`TableauHeun3`](@ref).

$(reference(Val(:Heun3)))
"""
struct Heun3 <: ERKMethod end

"""
Explicit Runge-Kutta method with [`TableauKutta3`](@ref).

$(reference(Val(:Kutta3)))
"""
struct Kutta3 <: ERKMethod end

"""
Explicit Runge-Kutta method with [`TableauRalston2`](@ref).

$(reference(Val(:Ralston2)))
"""
struct Ralston2 <: ERKMethod end

"""
Explicit Runge-Kutta method with [`TableauRalston3`](@ref).

$(reference(Val(:Ralston3)))
"""
struct Ralston3 <: ERKMethod end

"""
Explicit Runge-Kutta method with [`TableauRK4`](@ref).

$(reference(Val(:RK4)))
"""
struct RK4 <: ERKMethod end

"""
Explicit Runge-Kutta method with [`TableauRK21`](@ref).

$(reference(Val(:RK21)))
"""
struct RK21 <: ERKMethod end

"""
Explicit Runge-Kutta method with [`TableauRK22`](@ref).

$(reference(Val(:RK22)))
"""
struct RK22 <: ERKMethod end

"""
Explicit Runge-Kutta method with [`TableauRK31`](@ref).

$(reference(Val(:RK31)))
"""
struct RK31 <: ERKMethod end

"""
Explicit Runge-Kutta method with [`TableauRK32`](@ref).

$(reference(Val(:RK32)))
"""
struct RK32 <: ERKMethod end

"""
Explicit Runge-Kutta method with [`TableauRK41`](@ref).

$(reference(Val(:RK41)))
"""
struct RK41 <: ERKMethod end

"""
Explicit Runge-Kutta method with [`TableauRK42`](@ref).

$(reference(Val(:RK42)))
"""
struct RK42 <: ERKMethod end

"""
Explicit Runge-Kutta method with [`TableauRK416`](@ref).

$(reference(Val(:RK416)))
"""
struct RK416 <: ERKMethod end

"""
Explicit Runge-Kutta method with [`TableauRK438`](@ref).

$(reference(Val(:RK438)))
"""
struct RK438 <: ERKMethod end

"""
Explicit Runge-Kutta method with [`TableauRK5`](@ref).

$(reference(Val(:RK5)))
"""
struct RK5 <: ERKMethod end

"""
Explicit Runge-Kutta method with [`TableauRunge2`](@ref).

$(reference(Val(:Runge2)))
"""
struct Runge2 <: ERKMethod end

"""
Explicit Runge-Kutta method with [`TableauSSPRK2`](@ref).

$(reference(Val(:SSPRK2)))
"""
struct SSPRK2 <: ERKMethod end

"""
Explicit Runge-Kutta method with [`TableauSSPRK3`](@ref).

$(reference(Val(:SSPRK3)))
"""
struct SSPRK3 <: ERKMethod end

GeometricBase.tableau(::ForwardEuler) = TableauForwardEuler()
GeometricBase.tableau(::ExplicitEulerRK) = TableauExplicitEuler()
GeometricBase.tableau(::ExplicitMidpoint) = TableauExplicitMidpoint()
GeometricBase.tableau(::Heun2) = TableauHeun2()
GeometricBase.tableau(::Heun3) = TableauHeun3()
GeometricBase.tableau(::Kutta3) = TableauKutta3()
GeometricBase.tableau(::Ralston2) = TableauRalston2()
GeometricBase.tableau(::Ralston3) = TableauRalston3()
GeometricBase.tableau(::RK21) = TableauRK21()
GeometricBase.tableau(::RK22) = TableauRK22()
GeometricBase.tableau(::RK31) = TableauRK31()
GeometricBase.tableau(::RK32) = TableauRK32()
GeometricBase.tableau(::RK4) = TableauRK4()
GeometricBase.tableau(::RK41) = TableauRK41()
GeometricBase.tableau(::RK42) = TableauRK42()
GeometricBase.tableau(::RK416) = TableauRK416()
GeometricBase.tableau(::RK438) = TableauRK438()
GeometricBase.tableau(::RK5) = TableauRK5()
GeometricBase.tableau(::Runge2) = TableauRunge2()
GeometricBase.tableau(::SSPRK2) = TableauSSPRK2()
GeometricBase.tableau(::SSPRK3) = TableauSSPRK3()


# Diagonally Implicit Runge-Kutta Methods

"""
Diagonally implicit Runge-Kutta method with [`TableauCrankNicolson`](@ref).

$(reference(Val(:CrankNicolson)))
"""
struct CrankNicolson <: DIRKMethod end

"""
Diagonally implicit Runge-Kutta method with [`TableauCrouzeix`](@ref).

$(reference(Val(:Crouzeix)))
"""
struct Crouzeix <: DIRKMethod end

"""
Diagonally implicit Runge-Kutta method with [`TableauKraaijevangerSpijker`](@ref).

$(reference(Val(:KraaijevangerSpijker)))
"""
struct KraaijevangerSpijker <: DIRKMethod end

"""
Diagonally implicit Runge-Kutta method with [`TableauQinZhang`](@ref).

$(reference(Val(:QinZhang)))
"""
struct QinZhang <: DIRKMethod end

GeometricBase.tableau(::CrankNicolson) = TableauCrankNicolson()
GeometricBase.tableau(::Crouzeix) = TableauCrouzeix()
GeometricBase.tableau(::KraaijevangerSpijker) = TableauKraaijevangerSpijker()
GeometricBase.tableau(::QinZhang) = TableauQinZhang()


# Fully Implicit Runge-Kutta Methods

"""
Fully implicit Runge-Kutta method with [`TableauBackwardEuler`](@ref).

$(reference(Val(:BackwardEuler)))
"""
struct BackwardEuler <: IRKMethod end

"""
Fully implicit Runge-Kutta method with [`TableauImplicitEuler`](@ref).

$(reference(Val(:ImplicitEuler)))
"""
struct ImplicitEulerRK <: IRKMethod end

"""
Fully implicit Runge-Kutta method with [`TableauImplicitMidpoint`](@ref).

$(reference(Val(:ImplicitMidpoint)))
"""
struct ImplicitMidpoint <: IRKMethod end

"""
Fully implicit Runge-Kutta method with [`TableauSRK3`](@ref).

$(reference(Val(:SRK3)))
"""
struct SRK3 <: IRKMethod end

GeometricBase.tableau(::BackwardEuler) = TableauBackwardEuler()
GeometricBase.tableau(::ImplicitEulerRK) = TableauImplicitEuler()
GeometricBase.tableau(::ImplicitMidpoint) = TableauImplicitMidpoint()
GeometricBase.tableau(::SRK3) = TableauSRK3()

"""
Fully implicit Runge-Kutta method with [`TableauGauss`](@ref).

$(reference(Val(:Gauss)))
"""
struct Gauss <: IRKMethod
    s::Int
end

"""
Runge-Kutta method with [`TableauLobattoIII`](@ref).

$(reference(Val(:LobattoIII)))
"""
struct LobattoIII <: RKMethod
    s::Int
end

"""
Runge-Kutta method with [`TableauLobattoIIIA`](@ref).

$(reference(Val(:LobattoIIIA)))
"""
struct LobattoIIIA <: RKMethod
    s::Int
end

"""
Runge-Kutta method with [`TableauLobattoIIIB`](@ref).

$(reference(Val(:LobattoIIIB)))
"""
struct LobattoIIIB <: RKMethod
    s::Int
end

"""
Runge-Kutta method with [`TableauLobattoIIIC`](@ref).

$(reference(Val(:LobattoIIIC)))
"""
struct LobattoIIIC <: RKMethod
    s::Int
end

"""
Runge-Kutta method with [`TableauLobattoIIID`](@ref).

$(reference(Val(:LobattoIIID)))
"""
struct LobattoIIID <: RKMethod
    s::Int
end

"""
Runge-Kutta method with [`TableauLobattoIIIE`](@ref).

$(reference(Val(:LobattoIIIE)))
"""
struct LobattoIIIE <: RKMethod
    s::Int
end

"""
Runge-Kutta method with [`TableauLobattoIIIF`](@ref).

$(reference(Val(:LobattoIIIF)))
"""
struct LobattoIIIF <: RKMethod
    s::Int
end

"""
Runge-Kutta method with [`TableauLobattoIIIF̄`](@ref).
"""
struct LobattoIIIF̄ <: RKMethod
    s::Int
end

"""
Runge-Kutta method with [`TableauLobattoIIIG`](@ref).
"""
struct LobattoIIIG <: RKMethod
    s::Int
end

"""
Runge-Kutta method with [`TableauRadauIA`](@ref).

$(reference(Val(:RadauIA)))
"""
struct RadauIA <: RKMethod
    s::Int
end

"""
Runge-Kutta method with [`TableauRadauIB`](@ref).

$(reference(Val(:RadauIB)))
"""
struct RadauIB <: RKMethod
    s::Int
end

"""
Runge-Kutta method with [`TableauRadauIIA`](@ref).

$(reference(Val(:RadauIIA)))
"""
struct RadauIIA <: RKMethod
    s::Int
end

"""
Runge-Kutta method with [`TableauRadauIIB`](@ref).

$(reference(Val(:RadauIIB)))
"""
struct RadauIIB <: RKMethod
    s::Int
end

GeometricBase.tableau(method::Gauss) = TableauGauss(method.s)
GeometricBase.tableau(method::LobattoIII) = TableauLobattoIII(method.s)
GeometricBase.tableau(method::LobattoIIIA) = TableauLobattoIIIA(method.s)
GeometricBase.tableau(method::LobattoIIIB) = TableauLobattoIIIB(method.s)
GeometricBase.tableau(method::LobattoIIIC) = TableauLobattoIIIC(method.s)
GeometricBase.tableau(method::LobattoIIID) = TableauLobattoIIID(method.s)
GeometricBase.tableau(method::LobattoIIIE) = TableauLobattoIIIE(method.s)
GeometricBase.tableau(method::LobattoIIIF) = TableauLobattoIIIF(method.s)
GeometricBase.tableau(method::LobattoIIIF̄) = TableauLobattoIIIF̄(method.s)
GeometricBase.tableau(method::LobattoIIIG) = TableauLobattoIIIG(method.s)
GeometricBase.tableau(method::RadauIA) = TableauRadauIA(method.s)
GeometricBase.tableau(method::RadauIB) = TableauRadauIB(method.s)
GeometricBase.tableau(method::RadauIIA) = TableauRadauIIA(method.s)
GeometricBase.tableau(method::RadauIIB) = TableauRadauIIB(method.s)

GeometricBase.order(::Type{Gauss}) = "2s"
GeometricBase.order(::Type{LobattoIII}) = "2s-2"
GeometricBase.order(::Type{LobattoIIIA}) = "2s-2"
GeometricBase.order(::Type{LobattoIIIB}) = "2s-2"
GeometricBase.order(::Type{LobattoIIIC}) = "2s-2"
GeometricBase.order(::Type{LobattoIIID}) = "2s-2"
GeometricBase.order(::Type{LobattoIIIE}) = "2s-2"
GeometricBase.order(::Type{LobattoIIIF}) = "2s"
GeometricBase.order(::Type{LobattoIIIF̄}) = "2s"
GeometricBase.order(::Type{LobattoIIIG}) = "2s"
GeometricBase.order(::Type{RadauIA}) = "2s-1"
GeometricBase.order(::Type{RadauIB}) = "2s-1"
GeometricBase.order(::Type{RadauIIA}) = "2s-1"
GeometricBase.order(::Type{RadauIIB}) = "2s-1"

isexplicit(::Type{Gauss}) = false
isexplicit(::Type{LobattoIII}) = false
isexplicit(::Type{LobattoIIIA}) = false
isexplicit(::Type{LobattoIIIB}) = false
isexplicit(::Type{LobattoIIIC}) = false
isexplicit(::Type{LobattoIIID}) = false
isexplicit(::Type{LobattoIIIE}) = false
isexplicit(::Type{LobattoIIIF}) = false
isexplicit(::Type{LobattoIIIF̄}) = false
isexplicit(::Type{LobattoIIIG}) = false
isexplicit(::Type{RadauIA}) = false
isexplicit(::Type{RadauIB}) = false
isexplicit(::Type{RadauIIA}) = false
isexplicit(::Type{RadauIIB}) = false

isimplicit(::Type{Gauss}) = true
isimplicit(::Type{LobattoIII}) = true
isimplicit(::Type{LobattoIIIA}) = true
isimplicit(::Type{LobattoIIIB}) = true
isimplicit(::Type{LobattoIIIC}) = true
isimplicit(::Type{LobattoIIID}) = true
isimplicit(::Type{LobattoIIIE}) = true
isimplicit(::Type{LobattoIIIF}) = true
isimplicit(::Type{LobattoIIIF̄}) = true
isimplicit(::Type{LobattoIIIG}) = true
isimplicit(::Type{RadauIA}) = true
isimplicit(::Type{RadauIB}) = true
isimplicit(::Type{RadauIIA}) = true
isimplicit(::Type{RadauIIB}) = true

issymmetric(::Type{Gauss}) = true
issymmetric(::Type{LobattoIII}) = false
issymmetric(::Type{LobattoIIIA}) = true
issymmetric(::Type{LobattoIIIB}) = true
issymmetric(::Type{LobattoIIIC}) = false
issymmetric(::Type{LobattoIIID}) = true
issymmetric(::Type{LobattoIIIE}) = true
issymmetric(::Type{LobattoIIIF}) = true
issymmetric(::Type{LobattoIIIF̄}) = true
issymmetric(::Type{LobattoIIIG}) = true
issymmetric(::Type{RadauIA}) = false
issymmetric(::Type{RadauIB}) = false
issymmetric(::Type{RadauIIA}) = false
issymmetric(::Type{RadauIIB}) = false

issymplectic(::Type{Gauss}) = true
issymplectic(::Type{LobattoIII}) = false
issymplectic(::Type{LobattoIIIA}) = false
issymplectic(::Type{LobattoIIIB}) = false
issymplectic(::Type{LobattoIIIC}) = false
issymplectic(::Type{LobattoIIID}) = true
issymplectic(::Type{LobattoIIIE}) = true
issymplectic(::Type{LobattoIIIF}) = false
issymplectic(::Type{LobattoIIIF̄}) = false
issymplectic(::Type{LobattoIIIG}) = true
issymplectic(::Type{RadauIA}) = false
issymplectic(::Type{RadauIB}) = false
issymplectic(::Type{RadauIIA}) = false
issymplectic(::Type{RadauIIB}) = false


# Partitioned Runge-Kutta Methods


"""
    SymplecticEulerA

Symplectic Euler method using explicit Euler for q and implicit Euler for p.
"""
struct SymplecticEulerA <: EPRKMethod end

"""
    SymplecticEulerB

Symplectic Euler method using implicit Euler for q and explicit Euler for p.
"""
struct SymplecticEulerB <: EPRKMethod end
 
GeometricBase.tableau(::SymplecticEulerA) = PartitionedTableau(:SymplecticEulerA, TableauExplicitEuler(), TableauImplicitEuler())
GeometricBase.tableau(::SymplecticEulerB) = PartitionedTableau(:SymplecticEulerB, TableauImplicitEuler(), TableauExplicitEuler())

GeometricBase.order(::Type{SymplecticEulerA}) = 1
GeometricBase.order(::Type{SymplecticEulerB}) = 1

issymplectic(::Type{SymplecticEulerA}) = true
issymplectic(::Type{SymplecticEulerB}) = true


"""
Partitioned Runge-Kutta method [`TableauGauss`](@ref) for both ``q`` and ``p``.

$(reference(Val(:Gauss)))
"""
struct PartitionedGauss <: IPRKMethod
    s::Int
end


"Partitioned Runge-Kutta method with [`TableauLobattoIIIA`](@ref) for ``q`` and [`TableauLobattoIIIB`](@ref) for ``p``."
struct LobattoIIIAIIIB <: PRKMethod
    s::Int
end

"Partitioned Runge-Kutta method with [`TableauLobattoIIIB`](@ref) for ``q`` and [`TableauLobattoIIIA`](@ref) for ``p``."
struct LobattoIIIBIIIA <: PRKMethod
    s::Int
end

"Partitioned Runge-Kutta method with [`TableauLobattoIIIA`](@ref) for ``q`` and [`TableauLobattoIIIĀ`](@ref) for ``p``."
struct LobattoIIIAIIIĀ <: PRKMethod
    s::Int
end

"Partitioned Runge-Kutta method with [`TableauLobattoIIIB`](@ref) for ``q`` and [`TableauLobattoIIIB̄`](@ref) for ``p``."
struct LobattoIIIBIIIB̄ <: PRKMethod
    s::Int
end

"Partitioned Runge-Kutta method with [`TableauLobattoIIIC`](@ref) for ``q`` and [`TableauLobattoIIIC̄`](@ref) for ``p``."
struct LobattoIIICIIIC̄ <: PRKMethod
    s::Int
end

"Partitioned Runge-Kutta method with [`TableauLobattoIIIC̄`](@ref) for ``q`` and [`TableauLobattoIIIC`](@ref) for ``p``."
struct LobattoIIIC̄IIIC <: PRKMethod
    s::Int
end

"Partitioned Runge-Kutta method with [`TableauLobattoIIID`](@ref) for ``q`` and [`TableauLobattoIIID̄`](@ref) for ``p``."
struct LobattoIIIDIIID̄ <: PRKMethod
    s::Int
end

"Partitioned Runge-Kutta method with [`TableauLobattoIIIE`](@ref) for ``q`` and [`TableauLobattoIIIĒ`](@ref) for ``p``."
struct LobattoIIIEIIIĒ <: PRKMethod
    s::Int
end

"Partitioned Runge-Kutta method with [`TableauLobattoIIIF`](@ref) for ``q`` and [`TableauLobattoIIIF̄`](@ref) for ``p``."
struct LobattoIIIFIIIF̄ <: PRKMethod
    s::Int
end

"Partitioned Runge-Kutta method with [`TableauLobattoIIIF̄`](@ref) for ``q`` and [`TableauLobattoIIIF`](@ref) for ``p``."
struct LobattoIIIF̄IIIF <: PRKMethod
    s::Int
end

"Partitioned Runge-Kutta method with [`TableauLobattoIIIG`](@ref) for ``q`` and [`TableauLobattoIIIḠ`](@ref) for ``p``."
struct LobattoIIIGIIIḠ <: PRKMethod
    s::Int
end


GeometricBase.tableau(method::PartitionedGauss) = PartitionedTableau(TableauGauss(method.s))
GeometricBase.tableau(method::LobattoIIIAIIIB) = TableauLobattoIIIAIIIB(method.s)
GeometricBase.tableau(method::LobattoIIIBIIIA) = TableauLobattoIIIBIIIA(method.s)
GeometricBase.tableau(method::LobattoIIIAIIIĀ) = TableauLobattoIIIAIIIĀ(method.s)
GeometricBase.tableau(method::LobattoIIIBIIIB̄) = TableauLobattoIIIBIIIB̄(method.s)
GeometricBase.tableau(method::LobattoIIICIIIC̄) = TableauLobattoIIICIIIC̄(method.s)
GeometricBase.tableau(method::LobattoIIIC̄IIIC) = TableauLobattoIIIC̄IIIC(method.s)
GeometricBase.tableau(method::LobattoIIIDIIID̄) = TableauLobattoIIIDIIID̄(method.s)
GeometricBase.tableau(method::LobattoIIIEIIIĒ) = TableauLobattoIIIEIIIĒ(method.s)
GeometricBase.tableau(method::LobattoIIIFIIIF̄) = TableauLobattoIIIFIIIF̄(method.s)
GeometricBase.tableau(method::LobattoIIIF̄IIIF) = TableauLobattoIIIF̄IIIF(method.s)
GeometricBase.tableau(method::LobattoIIIGIIIḠ) = TableauLobattoIIIGIIIḠ(method.s)

GeometricBase.order(::Type{PartitionedGauss}) = "2s"
GeometricBase.order(::Type{LobattoIIIAIIIB}) = "2s-2"
GeometricBase.order(::Type{LobattoIIIBIIIA}) = "2s-2"
GeometricBase.order(::Type{LobattoIIIAIIIĀ}) = "2s-2"
GeometricBase.order(::Type{LobattoIIIBIIIB̄}) = "2s-2"
GeometricBase.order(::Type{LobattoIIICIIIC̄}) = "2s-2"
GeometricBase.order(::Type{LobattoIIIC̄IIIC}) = "2s-2"
GeometricBase.order(::Type{LobattoIIIDIIID̄}) = "2s-2"
GeometricBase.order(::Type{LobattoIIIEIIIĒ}) = "2s-2"
GeometricBase.order(::Type{LobattoIIIFIIIF̄}) = "2s"
GeometricBase.order(::Type{LobattoIIIF̄IIIF}) = "2s"
GeometricBase.order(::Type{LobattoIIIGIIIḠ}) = "2s"

issymplectic(::Type{PartitionedGauss}) = true
issymplectic(::Type{LobattoIIIAIIIB}) = true
issymplectic(::Type{LobattoIIIBIIIA}) = true
issymplectic(::Type{LobattoIIIAIIIĀ}) = true
issymplectic(::Type{LobattoIIIBIIIB̄}) = true
issymplectic(::Type{LobattoIIICIIIC̄}) = true
issymplectic(::Type{LobattoIIIC̄IIIC}) = true
issymplectic(::Type{LobattoIIIDIIID̄}) = true
issymplectic(::Type{LobattoIIIEIIIĒ}) = true
issymplectic(::Type{LobattoIIIFIIIF̄}) = true
issymplectic(::Type{LobattoIIIF̄IIIF}) = true
issymplectic(::Type{LobattoIIIGIIIḠ}) = true
