using GeometricEquations.Tests.HarmonicOscillator

using GeometricIntegrators.Methods: tableau

ode  = odeproblem()
pode = podeproblem()
hode = hodeproblem()
iode = iodeproblem()
lode = lodeproblem()


@testset "$(rpad("Runge-Kutta methods",80))" begin

    @test typeof(Integrator(ode, ForwardEuler())) <: IntegratorERK
    @test typeof(Integrator(ode, ExplicitEulerRK())) <: IntegratorERK
    @test typeof(Integrator(ode, ExplicitMidpoint())) <: IntegratorERK
    @test typeof(Integrator(ode, Heun2())) <: IntegratorERK
    @test typeof(Integrator(ode, Heun3())) <: IntegratorERK
    @test typeof(Integrator(ode, Kutta3())) <: IntegratorERK
    @test typeof(Integrator(ode, Ralston2())) <: IntegratorERK
    @test typeof(Integrator(ode, Ralston3())) <: IntegratorERK
    @test typeof(Integrator(ode, RK4())) <: IntegratorERK
    @test typeof(Integrator(ode, RK21())) <: IntegratorERK
    @test typeof(Integrator(ode, RK22())) <: IntegratorERK
    @test typeof(Integrator(ode, RK31())) <: IntegratorERK
    @test typeof(Integrator(ode, RK32())) <: IntegratorERK
    @test typeof(Integrator(ode, RK41())) <: IntegratorERK
    @test typeof(Integrator(ode, RK42())) <: IntegratorERK
    @test typeof(Integrator(ode, RK416())) <: IntegratorERK
    @test typeof(Integrator(ode, RK438())) <: IntegratorERK
    @test typeof(Integrator(ode, RK5())) <: IntegratorERK
    @test typeof(Integrator(ode, Runge2())) <: IntegratorERK
    @test typeof(Integrator(ode, SSPRK2())) <: IntegratorERK
    @test typeof(Integrator(ode, SSPRK3())) <: IntegratorERK

    @test typeof(Integrator(ode, CrankNicolson())) <: IntegratorDIRK
    @test typeof(Integrator(ode, Crouzeix())) <: IntegratorDIRK
    @test typeof(Integrator(ode, KraaijevangerSpijker())) <: IntegratorDIRK
    @test typeof(Integrator(ode, QinZhang())) <: IntegratorDIRK

    @test typeof(Integrator(ode, BackwardEuler())) <: IntegratorIRK
    @test typeof(Integrator(ode, ImplicitEulerRK())) <: IntegratorIRK
    @test typeof(Integrator(ode, ImplicitMidpoint())) <: IntegratorIRK
    @test typeof(Integrator(ode, SRK3())) <: IntegratorIRK

    @test typeof(Integrator(ode, Gauss(2))) <: IntegratorIRK
    @test typeof(Integrator(ode, LobattoIII(2))) <: IntegratorERK
    @test typeof(Integrator(ode, LobattoIIIA(2))) <: IntegratorDIRK
    @test typeof(Integrator(ode, LobattoIIIB(2))) <: IntegratorERK
    @test typeof(Integrator(ode, LobattoIIIC(2))) <: IntegratorIRK
    @test typeof(Integrator(ode, LobattoIIID(2))) <: IntegratorIRK
    @test typeof(Integrator(ode, LobattoIIIE(2))) <: IntegratorDIRK
    @test typeof(Integrator(ode, LobattoIIIF(2))) <: IntegratorIRK
    @test typeof(Integrator(ode, LobattoIIIF̄(2))) <: IntegratorIRK
    @test typeof(Integrator(ode, LobattoIIIG(2))) <: IntegratorIRK
    @test typeof(Integrator(ode, RadauIA(2))) <: IntegratorIRK
    @test typeof(Integrator(ode, RadauIB(2))) <: IntegratorIRK
    @test typeof(Integrator(ode, RadauIIA(2))) <: IntegratorIRK
    @test typeof(Integrator(ode, RadauIIB(2))) <: IntegratorIRK

    @test tableau(RK4()) == tableau(RK(TableauRK4()))

end


@testset "$(rpad("Partitioned Runge-Kutta methods",80))" begin

    @test typeof(Integrator(pode, SRK3())) <: IntegratorIPRK
    
    @test typeof(Integrator(pode, Gauss(2))) <: IntegratorIPRK
    @test typeof(Integrator(pode, RadauIA(2))) <: IntegratorIPRK
    @test typeof(Integrator(pode, RadauIB(2))) <: IntegratorIPRK
    @test typeof(Integrator(pode, RadauIIA(2))) <: IntegratorIPRK
    @test typeof(Integrator(pode, RadauIIB(2))) <: IntegratorIPRK

    @test typeof(Integrator(pode, LobattoIII(2))) <: IntegratorEPRK
    @test typeof(Integrator(pode, LobattoIIIA(2))) <: IntegratorIPRK
    @test typeof(Integrator(pode, LobattoIIIB(2))) <: IntegratorIPRK
    @test typeof(Integrator(pode, LobattoIIIC(2))) <: IntegratorIPRK
    @test typeof(Integrator(pode, LobattoIIID(2))) <: IntegratorIPRK
    @test typeof(Integrator(pode, LobattoIIIE(2))) <: IntegratorIPRK
    @test typeof(Integrator(pode, LobattoIIIF(2))) <: IntegratorIPRK
    @test typeof(Integrator(pode, LobattoIIIF̄(2))) <: IntegratorIPRK
    @test typeof(Integrator(pode, LobattoIIIG(2))) <: IntegratorIPRK
    
    @test typeof(Integrator(pode, LobattoIIIAIIIB(2))) <: IntegratorEPRK
    @test typeof(Integrator(pode, LobattoIIIBIIIA(2))) <: IntegratorEPRK
    @test typeof(Integrator(pode, LobattoIIIAIIIĀ(2))) <: IntegratorEPRK
    @test typeof(Integrator(pode, LobattoIIIBIIIB̄(2))) <: IntegratorEPRK
    @test typeof(Integrator(pode, LobattoIIICIIIC̄(2))) <: IntegratorIPRK
    @test typeof(Integrator(pode, LobattoIIIC̄IIIC(2))) <: IntegratorIPRK
    @test typeof(Integrator(pode, LobattoIIIDIIID̄(2))) <: IntegratorIPRK
    @test typeof(Integrator(pode, LobattoIIIEIIIĒ(2))) <: IntegratorIPRK
    @test typeof(Integrator(pode, LobattoIIIFIIIF̄(2))) <: IntegratorIPRK
    @test typeof(Integrator(pode, LobattoIIIF̄IIIF(2))) <: IntegratorIPRK
    @test typeof(Integrator(pode, LobattoIIIGIIIḠ(2))) <: IntegratorIPRK

end


@testset "$(rpad("Runge-Kutta methods for Implicit Equations",80))" begin

    @test typeof(Integrator(iode, ImplicitMidpoint())) <: IntegratorIRKimplicit
    @test typeof(Integrator(iode, SRK3())) <: IntegratorIRKimplicit
    @test typeof(Integrator(iode, Gauss(2))) <: IntegratorIRKimplicit

end

# @testset "$(rpad("Formal Lagrangian Runge-Kutta methods",80))" begin

#     @test typeof(Integrator(lode, FLRK(Gauss(1)))) <: IntegratorFLRK

# end
