using GeometricProblems.HarmonicOscillator
using RungeKutta.Tableaus

using GeometricIntegrators.Integrators: tableau

ode  = odeproblem()
pode = podeproblem()
hode = hodeproblem()
iode = iodeproblem()
lode = lodeproblem()


@testset "$(rpad("Runge-Kutta methods",80))" begin

    @test typeof(GeometricIntegrator(ode, ForwardEuler())) <: GeometricIntegrator{<:ERK}
    @test typeof(GeometricIntegrator(ode, ExplicitEulerRK())) <: GeometricIntegrator{<:ERK}
    @test typeof(GeometricIntegrator(ode, ExplicitMidpoint())) <: GeometricIntegrator{<:ERK}
    @test typeof(GeometricIntegrator(ode, Heun2())) <: GeometricIntegrator{<:ERK}
    @test typeof(GeometricIntegrator(ode, Heun3())) <: GeometricIntegrator{<:ERK}
    @test typeof(GeometricIntegrator(ode, Kutta3())) <: GeometricIntegrator{<:ERK}
    @test typeof(GeometricIntegrator(ode, Ralston2())) <: GeometricIntegrator{<:ERK}
    @test typeof(GeometricIntegrator(ode, Ralston3())) <: GeometricIntegrator{<:ERK}
    @test typeof(GeometricIntegrator(ode, RK4())) <: GeometricIntegrator{<:ERK}
    @test typeof(GeometricIntegrator(ode, RK21())) <: GeometricIntegrator{<:ERK}
    @test typeof(GeometricIntegrator(ode, RK22())) <: GeometricIntegrator{<:ERK}
    @test typeof(GeometricIntegrator(ode, RK31())) <: GeometricIntegrator{<:ERK}
    @test typeof(GeometricIntegrator(ode, RK32())) <: GeometricIntegrator{<:ERK}
    @test typeof(GeometricIntegrator(ode, RK41())) <: GeometricIntegrator{<:ERK}
    @test typeof(GeometricIntegrator(ode, RK42())) <: GeometricIntegrator{<:ERK}
    @test typeof(GeometricIntegrator(ode, RK416())) <: GeometricIntegrator{<:ERK}
    @test typeof(GeometricIntegrator(ode, RK438())) <: GeometricIntegrator{<:ERK}
    @test typeof(GeometricIntegrator(ode, RK5())) <: GeometricIntegrator{<:ERK}
    @test typeof(GeometricIntegrator(ode, Runge2())) <: GeometricIntegrator{<:ERK}
    @test typeof(GeometricIntegrator(ode, SSPRK2())) <: GeometricIntegrator{<:ERK}
    @test typeof(GeometricIntegrator(ode, SSPRK3())) <: GeometricIntegrator{<:ERK}

    @test typeof(GeometricIntegrator(ode, CrankNicolson())) <: GeometricIntegrator{<:DIRK}
    @test typeof(GeometricIntegrator(ode, Crouzeix())) <: GeometricIntegrator{<:DIRK}
    @test typeof(GeometricIntegrator(ode, KraaijevangerSpijker())) <: GeometricIntegrator{<:DIRK}
    @test typeof(GeometricIntegrator(ode, QinZhang())) <: GeometricIntegrator{<:DIRK}

    @test typeof(GeometricIntegrator(ode, BackwardEuler())) <: GeometricIntegrator{<:IRK}
    @test typeof(GeometricIntegrator(ode, ImplicitEulerRK())) <: GeometricIntegrator{<:IRK}
    @test typeof(GeometricIntegrator(ode, ImplicitMidpoint())) <: GeometricIntegrator{<:IRK}
    @test typeof(GeometricIntegrator(ode, SRK3())) <: GeometricIntegrator{<:IRK}

    @test typeof(GeometricIntegrator(ode, Gauss(2))) <: GeometricIntegrator{<:IRK}
    @test typeof(GeometricIntegrator(ode, LobattoIII(2))) <: GeometricIntegrator{<:ERK}
    @test typeof(GeometricIntegrator(ode, LobattoIIIA(2))) <: GeometricIntegrator{<:DIRK}
    @test typeof(GeometricIntegrator(ode, LobattoIIIB(2))) <: GeometricIntegrator{<:ERK}
    @test typeof(GeometricIntegrator(ode, LobattoIIIC(2))) <: GeometricIntegrator{<:IRK}
    @test typeof(GeometricIntegrator(ode, LobattoIIID(2))) <: GeometricIntegrator{<:IRK}
    @test typeof(GeometricIntegrator(ode, LobattoIIIE(2))) <: GeometricIntegrator{<:DIRK}
    @test typeof(GeometricIntegrator(ode, LobattoIIIF(2))) <: GeometricIntegrator{<:IRK}
    @test typeof(GeometricIntegrator(ode, LobattoIIIF̄(2))) <: GeometricIntegrator{<:IRK}
    @test typeof(GeometricIntegrator(ode, LobattoIIIG(2))) <: GeometricIntegrator{<:IRK}
    @test typeof(GeometricIntegrator(ode, RadauIA(2))) <: GeometricIntegrator{<:IRK}
    @test typeof(GeometricIntegrator(ode, RadauIB(2))) <: GeometricIntegrator{<:IRK}
    @test typeof(GeometricIntegrator(ode, RadauIIA(2))) <: GeometricIntegrator{<:IRK}
    @test typeof(GeometricIntegrator(ode, RadauIIB(2))) <: GeometricIntegrator{<:IRK}

    @test tableau(RK4()) == tableau(RK(TableauRK4()))

end


@testset "$(rpad("Partitioned Runge-Kutta methods",80))" begin

    @test typeof(GeometricIntegrator(pode, SRK3())) <: GeometricIntegrator{<:IPRK}
    
    @test typeof(GeometricIntegrator(pode, Gauss(2))) <: GeometricIntegrator{<:IPRK}
    @test typeof(GeometricIntegrator(pode, RadauIA(2))) <: GeometricIntegrator{<:IPRK}
    @test typeof(GeometricIntegrator(pode, RadauIB(2))) <: GeometricIntegrator{<:IPRK}
    @test typeof(GeometricIntegrator(pode, RadauIIA(2))) <: GeometricIntegrator{<:IPRK}
    @test typeof(GeometricIntegrator(pode, RadauIIB(2))) <: GeometricIntegrator{<:IPRK}

    @test typeof(GeometricIntegrator(pode, LobattoIII(2))) <: GeometricIntegrator{<:EPRK}
    @test typeof(GeometricIntegrator(pode, LobattoIIIA(2))) <: GeometricIntegrator{<:IPRK}
    @test typeof(GeometricIntegrator(pode, LobattoIIIB(2))) <: GeometricIntegrator{<:IPRK}
    @test typeof(GeometricIntegrator(pode, LobattoIIIC(2))) <: GeometricIntegrator{<:IPRK}
    @test typeof(GeometricIntegrator(pode, LobattoIIID(2))) <: GeometricIntegrator{<:IPRK}
    @test typeof(GeometricIntegrator(pode, LobattoIIIE(2))) <: GeometricIntegrator{<:IPRK}
    @test typeof(GeometricIntegrator(pode, LobattoIIIF(2))) <: GeometricIntegrator{<:IPRK}
    @test typeof(GeometricIntegrator(pode, LobattoIIIF̄(2))) <: GeometricIntegrator{<:IPRK}
    @test typeof(GeometricIntegrator(pode, LobattoIIIG(2))) <: GeometricIntegrator{<:IPRK}
    
    @test typeof(GeometricIntegrator(pode, LobattoIIIAIIIB(2))) <: GeometricIntegrator{<:EPRK}
    @test typeof(GeometricIntegrator(pode, LobattoIIIBIIIA(2))) <: GeometricIntegrator{<:EPRK}
    @test typeof(GeometricIntegrator(pode, LobattoIIIAIIIĀ(2))) <: GeometricIntegrator{<:EPRK}
    @test typeof(GeometricIntegrator(pode, LobattoIIIBIIIB̄(2))) <: GeometricIntegrator{<:EPRK}
    @test typeof(GeometricIntegrator(pode, LobattoIIICIIIC̄(2))) <: GeometricIntegrator{<:IPRK}
    @test typeof(GeometricIntegrator(pode, LobattoIIIC̄IIIC(2))) <: GeometricIntegrator{<:IPRK}
    @test typeof(GeometricIntegrator(pode, LobattoIIIDIIID̄(2))) <: GeometricIntegrator{<:IPRK}
    @test typeof(GeometricIntegrator(pode, LobattoIIIEIIIĒ(2))) <: GeometricIntegrator{<:IPRK}
    @test typeof(GeometricIntegrator(pode, LobattoIIIFIIIF̄(2))) <: GeometricIntegrator{<:IPRK}
    @test typeof(GeometricIntegrator(pode, LobattoIIIF̄IIIF(2))) <: GeometricIntegrator{<:IPRK}
    @test typeof(GeometricIntegrator(pode, LobattoIIIGIIIḠ(2))) <: GeometricIntegrator{<:IPRK}

end


@testset "$(rpad("Runge-Kutta methods for Implicit Equations",80))" begin

    @test typeof(GeometricIntegrator(iode, ImplicitMidpoint())) <: GeometricIntegrator{<:IRK}
    @test typeof(GeometricIntegrator(iode, SRK3())) <: GeometricIntegrator{<:IRK}
    @test typeof(GeometricIntegrator(iode, Gauss(2))) <: GeometricIntegrator{<:IRK}

end

# @testset "$(rpad("Formal Lagrangian Runge-Kutta methods",80))" begin

#     @test typeof(GeometricIntegrator(lode, FLRK(Gauss(1)))) <: IntegratorFLRK

# end
