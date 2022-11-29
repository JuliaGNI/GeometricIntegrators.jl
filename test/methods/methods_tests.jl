using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
using Test

using GeometricIntegrators.Methods: tableau

ode  = harmonic_oscillator_ode()
pode = harmonic_oscillator_pode()
hode = harmonic_oscillator_hode()
iode = harmonic_oscillator_iode()


@testset "$(rpad("Runge-Kutta methods",80))" begin

    @test typeof(Integrator(ode, ForwardEuler())) <: IntegratorERK
    @test typeof(Integrator(ode, ExplicitEuler())) <: IntegratorERK
    @test typeof(Integrator(ode, ExplicitMidpoint())) <: IntegratorERK
    @test typeof(Integrator(ode, Heun2())) <: IntegratorERK
    @test typeof(Integrator(ode, Heun3())) <: IntegratorERK
    @test typeof(Integrator(ode, Kutta3())) <: IntegratorERK
    @test typeof(Integrator(ode, Ralston2())) <: IntegratorERK
    @test typeof(Integrator(ode, Ralston3())) <: IntegratorERK
    @test typeof(Integrator(ode, RK4())) <: IntegratorERK
    @test typeof(Integrator(ode, RK416())) <: IntegratorERK
    @test typeof(Integrator(ode, RK438())) <: IntegratorERK
    @test typeof(Integrator(ode, Runge2())) <: IntegratorERK
    @test typeof(Integrator(ode, SSPRK2())) <: IntegratorERK
    @test typeof(Integrator(ode, SSPRK3())) <: IntegratorERK

    @test typeof(Integrator(ode, CrankNicolson())) <: IntegratorDIRK
    @test typeof(Integrator(ode, Crouzeix())) <: IntegratorDIRK
    @test typeof(Integrator(ode, KraaijevangerSpijker())) <: IntegratorDIRK
    @test typeof(Integrator(ode, QinZhang())) <: IntegratorDIRK

    @test typeof(Integrator(ode, BackwardEuler())) <: IntegratorFIRK
    @test typeof(Integrator(ode, ImplicitEuler())) <: IntegratorFIRK
    @test typeof(Integrator(ode, ImplicitMidpoint())) <: IntegratorFIRK
    @test typeof(Integrator(ode, SRK3())) <: IntegratorFIRK

    @test typeof(Integrator(ode, Gauss(2))) <: IntegratorFIRK
    @test typeof(Integrator(ode, LobattoIIIA(2))) <: IntegratorDIRK
    @test typeof(Integrator(ode, LobattoIIIB(2))) <: IntegratorERK
    @test typeof(Integrator(ode, LobattoIIIC(2))) <: IntegratorFIRK
    @test typeof(Integrator(ode, LobattoIIIC̄(2))) <: IntegratorERK
    @test typeof(Integrator(ode, LobattoIIID(2))) <: IntegratorFIRK
    @test typeof(Integrator(ode, LobattoIIIE(2))) <: IntegratorDIRK
    @test typeof(Integrator(ode, LobattoIIIF(2))) <: IntegratorFIRK
    @test typeof(Integrator(ode, LobattoIIIG(2))) <: IntegratorFIRK
    @test typeof(Integrator(ode, RadauIA(2))) <: IntegratorFIRK
    @test typeof(Integrator(ode, RadauIB(2))) <: IntegratorFIRK
    @test typeof(Integrator(ode, RadauIIA(2))) <: IntegratorFIRK
    @test typeof(Integrator(ode, RadauIIB(2))) <: IntegratorFIRK

    @test tableau(ForwardEuler()) == tableau(RK(TableauForwardEuler()))
    @test tableau(ExplicitEuler()) == tableau(RK(TableauExplicitEuler()))
    @test tableau(ExplicitMidpoint()) == tableau(RK(TableauExplicitMidpoint()))
    @test tableau(Heun2()) == tableau(RK(TableauHeun2()))
    @test tableau(Heun3()) == tableau(RK(TableauHeun3()))
    @test tableau(Kutta3()) == tableau(RK(TableauKutta3()))
    @test tableau(Ralston2()) == tableau(RK(TableauRalston2()))
    @test tableau(Ralston3()) == tableau(RK(TableauRalston3()))
    @test tableau(RK4()) == tableau(RK(TableauRK4()))
    @test tableau(RK416()) == tableau(RK(TableauRK416()))
    @test tableau(RK438()) == tableau(RK(TableauRK438()))
    @test tableau(Runge2()) == tableau(RK(TableauRunge2()))
    @test tableau(SSPRK2()) == tableau(RK(TableauSSPRK2()))
    @test tableau(SSPRK3()) == tableau(RK(TableauSSPRK3()))

end


@testset "$(rpad("Partitioned Runge-Kutta methods",80))" begin

    @test typeof(Integrator(pode, Gauss(2))) <: IntegratorIPRK
    @test typeof(Integrator(pode, SRK3())) <: IntegratorIPRK

    @test typeof(Integrator(pode, LobattoIIIAIIIB(2))) <: IntegratorEPRK
    @test typeof(Integrator(pode, LobattoIIIBIIIA(2))) <: IntegratorEPRK
    @test typeof(Integrator(pode, LobattoIIICIIIC̄(2))) <: IntegratorIPRK
    @test typeof(Integrator(pode, LobattoIIIC̄IIIC(2))) <: IntegratorIPRK

end


@testset "$(rpad("Runge-Kutta methods for Implicit Equations",80))" begin

    @test typeof(Integrator(iode, ImplicitMidpoint())) <: IntegratorFIRKimplicit
    @test typeof(Integrator(iode, SRK3())) <: IntegratorFIRKimplicit
    @test typeof(Integrator(iode, Gauss(2))) <: IntegratorFIRKimplicit

end
