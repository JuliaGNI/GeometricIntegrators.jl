
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


# Explicit Runge-Kutta Methods

struct ForwardEuler <: RKMethod end
struct ExplicitEuler <: RKMethod end
struct ExplicitMidpoint <: RKMethod end
struct Heun2 <: RKMethod end
struct Heun3 <: RKMethod end
struct Kutta3 <: RKMethod end
struct Ralston2 <: RKMethod end
struct Ralston3 <: RKMethod end
struct RK4 <: RKMethod end
struct RK416 <: RKMethod end
struct RK438 <: RKMethod end
struct Runge2 <: RKMethod end
struct SSPRK2 <: RKMethod end
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

struct CrankNicolson <: RKMethod end
struct Crouzeix <: RKMethod end
struct KraaijevangerSpijker <: RKMethod end
struct QinZhang <: RKMethod end

tableau(::CrankNicolson) = TableauCrankNicolson()
tableau(::Crouzeix) = TableauCrouzeix()
tableau(::KraaijevangerSpijker) = TableauKraaijevangerSpijker()
tableau(::QinZhang) = TableauQinZhang()


# Fully Implicit Runge-Kutta Methods

struct BackwardEuler <: RKMethod end
struct ImplicitEuler <: RKMethod end
struct ImplicitMidpoint <: RKMethod end
struct SRK3 <: RKMethod end

tableau(::BackwardEuler) = TableauBackwardEuler()
tableau(::ImplicitEuler) = TableauImplicitEuler()
tableau(::ImplicitMidpoint) = TableauImplicitMidpoint()
tableau(::SRK3) = TableauSRK3()


struct Gauss <: RKMethod
    s::Int
end

struct LobattoIIIA <: RKMethod
    s::Int
end

struct LobattoIIIB <: RKMethod
    s::Int
end

struct LobattoIIIC <: RKMethod
    s::Int
end

struct LobattoIIIC̄ <: RKMethod
    s::Int
end

struct LobattoIIID <: RKMethod
    s::Int
end

struct LobattoIIIE <: RKMethod
    s::Int
end

struct LobattoIIIF <: RKMethod
    s::Int
end

struct LobattoIIIG <: RKMethod
    s::Int
end

struct RadauIA <: RKMethod
    s::Int
end

struct RadauIB <: RKMethod
    s::Int
end

struct RadauIIA <: RKMethod
    s::Int
end

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
tableau(method::LobattoIIIG) = TableauLobattoIIIG(method.s)
tableau(method::RadauIA) = TableauRadauIA(method.s)
tableau(method::RadauIB) = TableauRadauIB(method.s)
tableau(method::RadauIIA) = TableauRadauIIA(method.s)
tableau(method::RadauIIB) = TableauRadauIIB(method.s)


# Partitioned Runge-Kutta Methods

struct LobattoIIIAIIIB <: PRKMethod
    s::Int
end

struct LobattoIIIBIIIA <: PRKMethod
    s::Int
end

struct LobattoIIICIIIC̄ <: PRKMethod
    s::Int
end

struct LobattoIIIC̄IIIC <: PRKMethod
    s::Int
end

tableau(method::LobattoIIIAIIIB) = TableauLobattoIIIAIIIB(method.s)
tableau(method::LobattoIIIBIIIA) = TableauLobattoIIIBIIIA(method.s)
tableau(method::LobattoIIICIIIC̄) = TableauLobattoIIICIIIC̄(method.s)
tableau(method::LobattoIIIC̄IIIC) = TableauLobattoIIIC̄IIIC(method.s)
