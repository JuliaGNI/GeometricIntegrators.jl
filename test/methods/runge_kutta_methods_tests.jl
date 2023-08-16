using GeometricProblems.HarmonicOscillator
using RungeKutta.Tableaus

using GeometricIntegrators.Methods: tableau

ode  = odeproblem()
pode = podeproblem()
hode = hodeproblem()
iode = iodeproblem()
lode = lodeproblem()


@testset "$(rpad("Runge-Kutta methods",80))" begin

    @test typeof(GeometricIntegrator(ode, ForwardEuler())) <: IntegratorERK
    @test typeof(GeometricIntegrator(ode, ExplicitEulerRK())) <: IntegratorERK
    @test typeof(GeometricIntegrator(ode, ExplicitMidpoint())) <: IntegratorERK
    @test typeof(GeometricIntegrator(ode, Heun2())) <: IntegratorERK
    @test typeof(GeometricIntegrator(ode, Heun3())) <: IntegratorERK
    @test typeof(GeometricIntegrator(ode, Kutta3())) <: IntegratorERK
    @test typeof(GeometricIntegrator(ode, Ralston2())) <: IntegratorERK
    @test typeof(GeometricIntegrator(ode, Ralston3())) <: IntegratorERK
    @test typeof(GeometricIntegrator(ode, RK4())) <: IntegratorERK
    @test typeof(GeometricIntegrator(ode, RK21())) <: IntegratorERK
    @test typeof(GeometricIntegrator(ode, RK22())) <: IntegratorERK
    @test typeof(GeometricIntegrator(ode, RK31())) <: IntegratorERK
    @test typeof(GeometricIntegrator(ode, RK32())) <: IntegratorERK
    @test typeof(GeometricIntegrator(ode, RK41())) <: IntegratorERK
    @test typeof(GeometricIntegrator(ode, RK42())) <: IntegratorERK
    @test typeof(GeometricIntegrator(ode, RK416())) <: IntegratorERK
    @test typeof(GeometricIntegrator(ode, RK438())) <: IntegratorERK
    @test typeof(GeometricIntegrator(ode, RK5())) <: IntegratorERK
    @test typeof(GeometricIntegrator(ode, Runge2())) <: IntegratorERK
    @test typeof(GeometricIntegrator(ode, SSPRK2())) <: IntegratorERK
    @test typeof(GeometricIntegrator(ode, SSPRK3())) <: IntegratorERK

    @test typeof(GeometricIntegrator(ode, CrankNicolson())) <: IntegratorDIRK
    @test typeof(GeometricIntegrator(ode, Crouzeix())) <: IntegratorDIRK
    @test typeof(GeometricIntegrator(ode, KraaijevangerSpijker())) <: IntegratorDIRK
    @test typeof(GeometricIntegrator(ode, QinZhang())) <: IntegratorDIRK

    @test typeof(GeometricIntegrator(ode, BackwardEuler())) <: IntegratorIRK
    @test typeof(GeometricIntegrator(ode, ImplicitEulerRK())) <: IntegratorIRK
    @test typeof(GeometricIntegrator(ode, ImplicitMidpoint())) <: IntegratorIRK
    @test typeof(GeometricIntegrator(ode, SRK3())) <: IntegratorIRK

    @test typeof(GeometricIntegrator(ode, Gauss(2))) <: IntegratorIRK
    @test typeof(GeometricIntegrator(ode, LobattoIII(2))) <: IntegratorERK
    @test typeof(GeometricIntegrator(ode, LobattoIIIA(2))) <: IntegratorDIRK
    @test typeof(GeometricIntegrator(ode, LobattoIIIB(2))) <: IntegratorERK
    @test typeof(GeometricIntegrator(ode, LobattoIIIC(2))) <: IntegratorIRK
    @test typeof(GeometricIntegrator(ode, LobattoIIID(2))) <: IntegratorIRK
    @test typeof(GeometricIntegrator(ode, LobattoIIIE(2))) <: IntegratorDIRK
    @test typeof(GeometricIntegrator(ode, LobattoIIIF(2))) <: IntegratorIRK
    @test typeof(GeometricIntegrator(ode, LobattoIIIF̄(2))) <: IntegratorIRK
    @test typeof(GeometricIntegrator(ode, LobattoIIIG(2))) <: IntegratorIRK
    @test typeof(GeometricIntegrator(ode, RadauIA(2))) <: IntegratorIRK
    @test typeof(GeometricIntegrator(ode, RadauIB(2))) <: IntegratorIRK
    @test typeof(GeometricIntegrator(ode, RadauIIA(2))) <: IntegratorIRK
    @test typeof(GeometricIntegrator(ode, RadauIIB(2))) <: IntegratorIRK

    @test tableau(RK4()) == tableau(RK(TableauRK4()))

end


@testset "$(rpad("Partitioned Runge-Kutta methods",80))" begin

    @test typeof(GeometricIntegrator(pode, SRK3())) <: IntegratorIPRK
    
    @test typeof(GeometricIntegrator(pode, Gauss(2))) <: IntegratorIPRK
    @test typeof(GeometricIntegrator(pode, RadauIA(2))) <: IntegratorIPRK
    @test typeof(GeometricIntegrator(pode, RadauIB(2))) <: IntegratorIPRK
    @test typeof(GeometricIntegrator(pode, RadauIIA(2))) <: IntegratorIPRK
    @test typeof(GeometricIntegrator(pode, RadauIIB(2))) <: IntegratorIPRK

    @test typeof(GeometricIntegrator(pode, LobattoIII(2))) <: IntegratorEPRK
    @test typeof(GeometricIntegrator(pode, LobattoIIIA(2))) <: IntegratorIPRK
    @test typeof(GeometricIntegrator(pode, LobattoIIIB(2))) <: IntegratorIPRK
    @test typeof(GeometricIntegrator(pode, LobattoIIIC(2))) <: IntegratorIPRK
    @test typeof(GeometricIntegrator(pode, LobattoIIID(2))) <: IntegratorIPRK
    @test typeof(GeometricIntegrator(pode, LobattoIIIE(2))) <: IntegratorIPRK
    @test typeof(GeometricIntegrator(pode, LobattoIIIF(2))) <: IntegratorIPRK
    @test typeof(GeometricIntegrator(pode, LobattoIIIF̄(2))) <: IntegratorIPRK
    @test typeof(GeometricIntegrator(pode, LobattoIIIG(2))) <: IntegratorIPRK
    
    @test typeof(GeometricIntegrator(pode, LobattoIIIAIIIB(2))) <: IntegratorEPRK
    @test typeof(GeometricIntegrator(pode, LobattoIIIBIIIA(2))) <: IntegratorEPRK
    @test typeof(GeometricIntegrator(pode, LobattoIIIAIIIĀ(2))) <: IntegratorEPRK
    @test typeof(GeometricIntegrator(pode, LobattoIIIBIIIB̄(2))) <: IntegratorEPRK
    @test typeof(GeometricIntegrator(pode, LobattoIIICIIIC̄(2))) <: IntegratorIPRK
    @test typeof(GeometricIntegrator(pode, LobattoIIIC̄IIIC(2))) <: IntegratorIPRK
    @test typeof(GeometricIntegrator(pode, LobattoIIIDIIID̄(2))) <: IntegratorIPRK
    @test typeof(GeometricIntegrator(pode, LobattoIIIEIIIĒ(2))) <: IntegratorIPRK
    @test typeof(GeometricIntegrator(pode, LobattoIIIFIIIF̄(2))) <: IntegratorIPRK
    @test typeof(GeometricIntegrator(pode, LobattoIIIF̄IIIF(2))) <: IntegratorIPRK
    @test typeof(GeometricIntegrator(pode, LobattoIIIGIIIḠ(2))) <: IntegratorIPRK

end


@testset "$(rpad("Runge-Kutta methods for Implicit Equations",80))" begin

    @test typeof(GeometricIntegrator(iode, ImplicitMidpoint())) <: IntegratorIRKimplicit
    @test typeof(GeometricIntegrator(iode, SRK3())) <: IntegratorIRKimplicit
    @test typeof(GeometricIntegrator(iode, Gauss(2))) <: IntegratorIRKimplicit

end

# @testset "$(rpad("Formal Lagrangian Runge-Kutta methods",80))" begin

#     @test typeof(GeometricIntegrator(lode, FLRK(Gauss(1)))) <: IntegratorFLRK

# end
