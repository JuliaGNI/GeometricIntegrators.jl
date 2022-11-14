
abstract type RKODEMethod <: ODEMethod end
abstract type RKPODEMethod <: PODEMethod end

tableau(method::RKODEMethod) = error("No tableau for Runge-Kutta method $(typeof(method)) provided")
tableau(method::RKPODEMethod) = error("No tableau for partitioned Runge-Kutta method $(typeof(method)) provided")

Integrators.Integrator(problem::ODEProblem, method::RKODEMethod; kwargs...) = Integrator(problem, tableau(method); kwargs...)
Integrators.Integrator(problem::Union{PODEProblem,HODEProblem}, method::RKODEMethod; kwargs...) = Integrator(problem, PartitionedTableau(tableau(method)); kwargs...)
Integrators.Integrator(problem::Union{PODEProblem,HODEProblem}, method::RKPODEMethod; kwargs...) = Integrator(problem, tableau(method); kwargs...)


# Explicit Runge-Kutta Methods

struct ForwardEuler <: RKODEMethod end
struct ExplicitEuler <: RKODEMethod end
struct ExplicitMidpoint <: RKODEMethod end
struct Heun2 <: RKODEMethod end
struct Heun3 <: RKODEMethod end
struct Kutta3 <: RKODEMethod end
struct Ralston2 <: RKODEMethod end
struct Ralston3 <: RKODEMethod end
struct RK4 <: RKODEMethod end
struct RK416 <: RKODEMethod end
struct RK438 <: RKODEMethod end
struct Runge2 <: RKODEMethod end
struct SSPRK2 <: RKODEMethod end
struct SSPRK3 <: RKODEMethod end

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

struct CrankNicolson <: RKODEMethod end
struct Crouzeix <: RKODEMethod end
struct KraaijevangerSpijker <: RKODEMethod end
struct QinZhang <: RKODEMethod end

tableau(::CrankNicolson) = TableauCrankNicolson()
tableau(::Crouzeix) = TableauCrouzeix()
tableau(::KraaijevangerSpijker) = TableauKraaijevangerSpijker()
tableau(::QinZhang) = TableauQinZhang()


# Fully Implicit Runge-Kutta Methods

struct BackwardEuler <: RKODEMethod end
struct ImplicitEuler <: RKODEMethod end
struct ImplicitMidpoint <: RKODEMethod end
struct SRK3 <: RKODEMethod end

tableau(::BackwardEuler) = TableauBackwardEuler()
tableau(::ImplicitEuler) = TableauImplicitEuler()
tableau(::ImplicitMidpoint) = TableauImplicitMidpoint()
tableau(::SRK3) = TableauSRK3()


struct Gauss <: RKODEMethod
    s::Int
end

struct LobattoIIIA <: RKODEMethod
    s::Int
end

struct LobattoIIIB <: RKODEMethod
    s::Int
end

struct LobattoIIIC <: RKODEMethod
    s::Int
end

struct LobattoIIIC̄ <: RKODEMethod
    s::Int
end

struct LobattoIIID <: RKODEMethod
    s::Int
end

struct LobattoIIIE <: RKODEMethod
    s::Int
end

struct LobattoIIIF <: RKODEMethod
    s::Int
end

struct LobattoIIIG <: RKODEMethod
    s::Int
end

struct RadauIA <: RKODEMethod
    s::Int
end

struct RadauIB <: RKODEMethod
    s::Int
end

struct RadauIIA <: RKODEMethod
    s::Int
end

struct RadauIIB <: RKODEMethod
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
tableau(method::LobattoIIIG) = TableauLobattoIIIG(method.s)
tableau(method::RadauIA) = TableauRadauIA(method.s)
tableau(method::RadauIB) = TableauRadauIB(method.s)
tableau(method::RadauIIA) = TableauRadauIIA(method.s)
tableau(method::RadauIIB) = TableauRadauIIB(method.s)


# Partitioned Runge-Kutta Methods

struct LobattoIIIAIIIB <: RKPODEMethod
    s::Int
end

struct LobattoIIIBIIIA <: RKPODEMethod
    s::Int
end

struct LobattoIIICIIIC̄ <: RKPODEMethod
    s::Int
end

struct LobattoIIIC̄IIIC <: RKPODEMethod
    s::Int
end

tableau(method::LobattoIIIAIIIB) = TableauLobattoIIIAIIIB(method.s)
tableau(method::LobattoIIIBIIIA) = TableauLobattoIIIBIIIA(method.s)
tableau(method::LobattoIIICIIIC̄) = TableauLobattoIIICIIIC̄(method.s)
tableau(method::LobattoIIIC̄IIIC) = TableauLobattoIIIC̄IIIC(method.s)
