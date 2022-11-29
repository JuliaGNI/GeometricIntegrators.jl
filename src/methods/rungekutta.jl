
abstract type RKMethod <: ODEMethod end
abstract type PRKMethod <: PODEMethod end

tableau(method::RKMethod) = error("No tableau for Runge-Kutta method $(typeof(method)) provided")
tableau(method::PRKMethod) = error("No tableau for partitioned Runge-Kutta method $(typeof(method)) provided")

ispodemethod(::RKMethod) = true
ishodemethod(::RKMethod) = true
isiodemethod(::RKMethod) = true
islodemethod(::RKMethod) = true
isiodemethod(::PRKMethod) = true
islodemethod(::PRKMethod) = true

Integrators.Integrator(problem::ODEProblem, method::RKMethod; kwargs...) = Integrator(problem, tableau(method); kwargs...)
Integrators.Integrator(problem::Union{PODEProblem,HODEProblem}, method::RKMethod; kwargs...) = Integrator(problem, PartitionedTableau(tableau(method)); kwargs...)
Integrators.Integrator(problem::Union{PODEProblem,HODEProblem}, method::PRKMethod; kwargs...) = Integrator(problem, tableau(method); kwargs...)
Integrators.Integrator(problem::Union{IODEProblem,LODEProblem}, method::RKMethod; kwargs...) = Integrator(problem, tableau(method); kwargs...)
Integrators.Integrator(problem::Union{IODEProblem,LODEProblem}, method::PRKMethod; kwargs...) = Integrator(problem, tableau(method); kwargs...)


# General Runge-Kutta Method

struct RK{TT} <: RKMethod
    tableau::TT

    function RK(tableau::TT) where {TT <: Tableau}
        new{TT}(tableau)
    end
end

tableau(method::RK) = method.tableau

Integrators.Integrator(problem::ODEProblem, method::RK; kwargs...) = Integrator(problem, tableau(method); kwargs...)


# Explicit Runge-Kutta Methods

"Explicit Runge-Kutta method with [`TableauForwardEuler`](@ref)."
struct ForwardEuler <: RKMethod end

"Explicit Runge-Kutta method with [`TableauExplicitEuler`](@ref)."
struct ExplicitEuler <: RKMethod end

"Explicit Runge-Kutta method with [`TableauExplicitMidpoint`](@ref)."
struct ExplicitMidpoint <: RKMethod end

"Explicit Runge-Kutta method with [`TableauHeun2`](@ref)."
struct Heun2 <: RKMethod end

"Explicit Runge-Kutta method with [`TableauHeun3`](@ref)."
struct Heun3 <: RKMethod end

"Explicit Runge-Kutta method with [`TableauKutta3`](@ref)."
struct Kutta3 <: RKMethod end

"Explicit Runge-Kutta method with [`TableauRalston2`](@ref)."
struct Ralston2 <: RKMethod end

"Explicit Runge-Kutta method with [`TableauRalston3`](@ref)."
struct Ralston3 <: RKMethod end

"Explicit Runge-Kutta method with [`TableauRK4`](@ref)."
struct RK4 <: RKMethod end

"Explicit Runge-Kutta method with [`TableauRK416`](@ref)."
struct RK416 <: RKMethod end

"Explicit Runge-Kutta method with [`TableauRK438`](@ref)."
struct RK438 <: RKMethod end

"Explicit Runge-Kutta method with [`TableauRunge2`](@ref)."
struct Runge2 <: RKMethod end

"Explicit Runge-Kutta method with [`TableauSSPRK2`](@ref)."
struct SSPRK2 <: RKMethod end

"Explicit Runge-Kutta method with [`TableauSSPRK3`](@ref)."
struct SSPRK3 <: RKMethod end

tableau(::ForwardEuler) = TableauForwardEuler()
tableau(::ExplicitEuler) = TableauExplicitEuler()
tableau(::ExplicitMidpoint) = TableauExplicitMidpoint()
tableau(::Heun2) = TableauHeun2()
tableau(::Heun3) = TableauHeun3()
tableau(::Kutta3) = TableauKutta3()
tableau(::Ralston2) = TableauRalston2()
tableau(::Ralston3) = TableauRalston3()
tableau(::RK4) = TableauRK4()
tableau(::RK416) = TableauRK416()
tableau(::RK438) = TableauRK438()
tableau(::Runge2) = TableauRunge2()
tableau(::SSPRK2) = TableauSSPRK2()
tableau(::SSPRK3) = TableauSSPRK3()


# Diagonally Implicit Runge-Kutta Methods

"Diagonally implicit Runge-Kutta method with [`TableauCrankNicolson`](@ref)."
struct CrankNicolson <: RKMethod end

"Diagonally implicit Runge-Kutta method with [`TableauCrouzeix`](@ref)."
struct Crouzeix <: RKMethod end

"Diagonally implicit Runge-Kutta method with [`TableauKraaijevangerSpijker`](@ref)."
struct KraaijevangerSpijker <: RKMethod end

"Diagonally implicit Runge-Kutta method with [`TableauQinZhang`](@ref)."
struct QinZhang <: RKMethod end

tableau(::CrankNicolson) = TableauCrankNicolson()
tableau(::Crouzeix) = TableauCrouzeix()
tableau(::KraaijevangerSpijker) = TableauKraaijevangerSpijker()
tableau(::QinZhang) = TableauQinZhang()


# Fully Implicit Runge-Kutta Methods

"Fully implicit Runge-Kutta method with [`TableauBackwardEuler`](@ref)."
struct BackwardEuler <: RKMethod end

"Fully implicit Runge-Kutta method with [`TableauImplicitEuler`](@ref)."
struct ImplicitEuler <: RKMethod end

"Fully implicit Runge-Kutta method with [`TableauImplicitMidpoint`](@ref)."
struct ImplicitMidpoint <: RKMethod end

"Fully implicit Runge-Kutta method with [`TableauSRK3`](@ref)."
struct SRK3 <: RKMethod end

tableau(::BackwardEuler) = TableauBackwardEuler()
tableau(::ImplicitEuler) = TableauImplicitEuler()
tableau(::ImplicitMidpoint) = TableauImplicitMidpoint()
tableau(::SRK3) = TableauSRK3()

"Fully implicit Runge-Kutta method with [`TableauGauss`](@ref)."
struct Gauss <: RKMethod
    s::Int
end

"Runge-Kutta method with [`TableauLobattoIIIA`](@ref)."
struct LobattoIIIA <: RKMethod
    s::Int
end

"Runge-Kutta method with [`TableauLobattoIIIB`](@ref)."
struct LobattoIIIB <: RKMethod
    s::Int
end

"Runge-Kutta method with [`TableauLobattoIIIC`](@ref)."
struct LobattoIIIC <: RKMethod
    s::Int
end

"Runge-Kutta method with [`TableauLobattoIIIC̄`](@ref)."
struct LobattoIIIC̄ <: RKMethod
    s::Int
end

"Runge-Kutta method with [`TableauLobattoIIID`](@ref)."
struct LobattoIIID <: RKMethod
    s::Int
end

"Runge-Kutta method with [`TableauLobattoIIIE`](@ref)."
struct LobattoIIIE <: RKMethod
    s::Int
end

"Runge-Kutta method with [`TableauLobattoIIIF`](@ref)."
struct LobattoIIIF <: RKMethod
    s::Int
end

"Runge-Kutta method with [`TableauLobattoIIIF̄`](@ref)."
struct LobattoIIIF̄ <: RKMethod
    s::Int
end

"Runge-Kutta method with [`TableauLobattoIIIG`](@ref)."
struct LobattoIIIG <: RKMethod
    s::Int
end

"Runge-Kutta method with [`TableauRadauIA`](@ref)."
struct RadauIA <: RKMethod
    s::Int
end

"Runge-Kutta method with [`TableauRadauIB`](@ref)."
struct RadauIB <: RKMethod
    s::Int
end

"Runge-Kutta method with [`TableauRadauIIA`](@ref)."
struct RadauIIA <: RKMethod
    s::Int
end

"Runge-Kutta method with [`TableauRadauIIB`](@ref)."
struct RadauIIB <: RKMethod
    s::Int
end

tableau(method::Gauss) = TableauGauss(method.s)
tableau(method::LobattoIIIA) = TableauLobattoIIIA(method.s)
tableau(method::LobattoIIIB) = TableauLobattoIIIB(method.s)
tableau(method::LobattoIIIC) = TableauLobattoIIIC(method.s)
tableau(method::LobattoIIIC̄) = TableauLobattoIIIC̄(method.s)
tableau(method::LobattoIIID) = TableauLobattoIIID(method.s)
tableau(method::LobattoIIIE) = TableauLobattoIIIE(method.s)
tableau(method::LobattoIIIF) = TableauLobattoIIIF(method.s)
tableau(method::LobattoIIIF̄) = TableauLobattoIIIF̄(method.s)
tableau(method::LobattoIIIG) = TableauLobattoIIIG(method.s)
tableau(method::RadauIA) = TableauRadauIA(method.s)
tableau(method::RadauIB) = TableauRadauIB(method.s)
tableau(method::RadauIIA) = TableauRadauIIA(method.s)
tableau(method::RadauIIB) = TableauRadauIIB(method.s)


# Partitioned Runge-Kutta Methods

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


tableau(method::LobattoIIIAIIIB) = TableauLobattoIIIAIIIB(method.s)
tableau(method::LobattoIIIBIIIA) = TableauLobattoIIIBIIIA(method.s)
tableau(method::LobattoIIIAIIIĀ) = TableauLobattoIIIAIIIĀ(method.s)
tableau(method::LobattoIIIBIIIB̄) = TableauLobattoIIIBIIIB̄(method.s)
tableau(method::LobattoIIICIIIC̄) = TableauLobattoIIICIIIC̄(method.s)
tableau(method::LobattoIIIC̄IIIC) = TableauLobattoIIIC̄IIIC(method.s)
tableau(method::LobattoIIIDIIID̄) = TableauLobattoIIIDIIID̄(method.s)
tableau(method::LobattoIIIEIIIĒ) = TableauLobattoIIIEIIIĒ(method.s)
tableau(method::LobattoIIIFIIIF̄) = TableauLobattoIIIFIIIF̄(method.s)
tableau(method::LobattoIIIF̄IIIF) = TableauLobattoIIIF̄IIIF(method.s)
tableau(method::LobattoIIIGIIIḠ) = TableauLobattoIIIGIIIḠ(method.s)