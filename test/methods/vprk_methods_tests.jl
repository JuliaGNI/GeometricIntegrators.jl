using GeometricIntegrators
using GeometricIntegrators.Integrators.VPRK
using GeometricProblems.HarmonicOscillator
using Test


iode = harmonic_oscillator_iode()
# lode = harmonic_oscillator_lode()


@testset "$(rpad("Variational Partitioned Runge-Kutta Methods",80))" begin

    @test typeof(VPSRK3()) <: VPRKMethod

    @test typeof(VPRKGauss(1)) <: VPRKMethod
    @test typeof(VPRKGauss(2)) <: VPRKMethod
    @test typeof(VPRKGauss(3)) <: VPRKMethod
    @test typeof(VPRKGauss(4)) <: VPRKMethod

    @test typeof(VPRKRadauIIA(2)) <: VPRKMethod
    @test typeof(VPRKRadauIIA(3)) <: VPRKMethod
    @test typeof(VPRKRadauIIA(4)) <: VPRKMethod
    @test typeof(VPRKRadauIIA(5)) <: VPRKMethod

    @test typeof(VPRKRadauIIB(2)) <: VPRKMethod
    @test typeof(VPRKRadauIIB(3)) <: VPRKMethod
    @test typeof(VPRKRadauIIB(4)) <: VPRKMethod
    @test typeof(VPRKRadauIIB(5)) <: VPRKMethod

    @test typeof(VPRKLobattoIIIA(2)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIA(3)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIA(4)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIA(5)) <: VPRKMethod

    @test typeof(VPRKLobattoIIIB(2)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIB(3)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIB(4)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIB(5)) <: VPRKMethod

    @test typeof(VPRKLobattoIIIC(2)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIC(3)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIC(4)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIC(5)) <: VPRKMethod

    @test typeof(VPRKLobattoIIIC̄(2)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIC̄(3)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIC̄(4)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIC̄(5)) <: VPRKMethod

    @test typeof(VPRKLobattoIIID(2)) <: VPRKMethod
    @test typeof(VPRKLobattoIIID(3)) <: VPRKMethod
    @test typeof(VPRKLobattoIIID(4)) <: VPRKMethod
    @test typeof(VPRKLobattoIIID(5)) <: VPRKMethod

    @test typeof(VPRKLobattoIIIE(2)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIE(3)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIE(4)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIE(5)) <: VPRKMethod

    @test typeof(VPRKLobattoIIIF(2)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIF(3)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIF(4)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIF(5)) <: VPRKMethod

    @test typeof(VPRKLobattoIIIG(2)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIG(3)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIG(4)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIG(5)) <: VPRKMethod

    @test typeof(VPRKLobattoIIIAIIIB(2)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIAIIIB(3)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIAIIIB(4)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIAIIIB(5)) <: VPRKMethod

    @test typeof(VPRKLobattoIIIBIIIA(2)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIBIIIA(3)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIBIIIA(4)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIBIIIA(5)) <: VPRKMethod

    @test typeof(VPRKLobattoIIIAIIIĀ(2)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIAIIIĀ(3)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIAIIIĀ(4)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIAIIIĀ(5)) <: VPRKMethod

    @test typeof(VPRKLobattoIIIBIIIB̄(2)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIBIIIB̄(3)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIBIIIB̄(4)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIBIIIB̄(5)) <: VPRKMethod

    @test typeof(VPRKLobattoIIICIIIC̄(2)) <: VPRKMethod
    @test typeof(VPRKLobattoIIICIIIC̄(3)) <: VPRKMethod
    @test typeof(VPRKLobattoIIICIIIC̄(4)) <: VPRKMethod
    @test typeof(VPRKLobattoIIICIIIC̄(5)) <: VPRKMethod

    @test typeof(VPRKLobattoIIIC̄IIIC(2)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIC̄IIIC(3)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIC̄IIIC(4)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIC̄IIIC(5)) <: VPRKMethod

    @test typeof(VPRKLobattoIIIDIIID̄(2)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIDIIID̄(3)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIDIIID̄(4)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIDIIID̄(5)) <: VPRKMethod

    @test typeof(VPRKLobattoIIIEIIIĒ(2)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIEIIIĒ(3)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIEIIIĒ(4)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIEIIIĒ(5)) <: VPRKMethod

    @test typeof(VPRKLobattoIIIFIIIF̄(2)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIFIIIF̄(3)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIFIIIF̄(4)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIFIIIF̄(5)) <: VPRKMethod

    @test typeof(VPRKLobattoIIIF̄IIIF(2)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIF̄IIIF(3)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIF̄IIIF(4)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIF̄IIIF(5)) <: VPRKMethod

    @test typeof(VPRKLobattoIIIGIIIḠ(2)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIGIIIḠ(3)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIGIIIḠ(4)) <: VPRKMethod
    @test typeof(VPRKLobattoIIIGIIIḠ(5)) <: VPRKMethod


    @test VPSRK3()           == VPSRK3(Float64)
    @test VPRKGauss(1)       == VPRKGauss(Float64, 1)
    @test VPRKRadauIIA(2)    == VPRKRadauIIA(Float64, 2)
    @test VPRKRadauIIB(2)    == VPRKRadauIIB(Float64, 2)
    @test VPRKLobattoIIIA(2) == VPRKLobattoIIIA(Float64, 2)
    @test VPRKLobattoIIIB(2) == VPRKLobattoIIIB(Float64, 2)
    @test VPRKLobattoIIIC(2) == VPRKLobattoIIIC(Float64, 2)
    @test VPRKLobattoIIIC̄(2) == VPRKLobattoIIIC̄(Float64, 2)
    @test VPRKLobattoIIID(2) == VPRKLobattoIIID(Float64, 2)
    @test VPRKLobattoIIIE(2) == VPRKLobattoIIIE(Float64, 2)
    @test VPRKLobattoIIIF(2) == VPRKLobattoIIIF(Float64, 2)
    @test VPRKLobattoIIIG(2) == VPRKLobattoIIIG(Float64, 2)


    @test typeof(Integrator(iode, VPSRK3())) <: IntegratorVPRK
    @test typeof(Integrator(iode, VPRKGauss(2))) <: IntegratorVPRK

    @test typeof(Integrator(iode, VPRKLobattoIIIA(2))) <: IntegratorVPRK
    @test typeof(Integrator(iode, VPRKLobattoIIIB(2))) <: IntegratorVPRK
    @test typeof(Integrator(iode, VPRKLobattoIIIC(2))) <: IntegratorVPRK
    @test typeof(Integrator(iode, VPRKLobattoIIIC̄(2))) <: IntegratorVPRK
    @test typeof(Integrator(iode, VPRKLobattoIIID(2))) <: IntegratorVPRK
    @test typeof(Integrator(iode, VPRKLobattoIIIE(2))) <: IntegratorVPRK
    @test typeof(Integrator(iode, VPRKLobattoIIIF(2))) <: IntegratorVPRK
    @test typeof(Integrator(iode, VPRKLobattoIIIG(2))) <: IntegratorVPRK

    @test typeof(Integrator(iode, VPRKLobattoIIIAIIIB(2))) <: IntegratorVPRK
    @test typeof(Integrator(iode, VPRKLobattoIIIBIIIA(2))) <: IntegratorVPRK
    @test typeof(Integrator(iode, VPRKLobattoIIIAIIIĀ(2))) <: IntegratorVPRK
    @test typeof(Integrator(iode, VPRKLobattoIIIBIIIB̄(2))) <: IntegratorVPRK
    @test typeof(Integrator(iode, VPRKLobattoIIICIIIC̄(2))) <: IntegratorVPRK
    @test typeof(Integrator(iode, VPRKLobattoIIIC̄IIIC(2))) <: IntegratorVPRK
    @test typeof(Integrator(iode, VPRKLobattoIIIDIIID̄(2))) <: IntegratorVPRK
    @test typeof(Integrator(iode, VPRKLobattoIIIEIIIĒ(2))) <: IntegratorVPRK
    @test typeof(Integrator(iode, VPRKLobattoIIIFIIIF̄(2))) <: IntegratorVPRK
    @test typeof(Integrator(iode, VPRKLobattoIIIF̄IIIF(2))) <: IntegratorVPRK
    @test typeof(Integrator(iode, VPRKLobattoIIIGIIIḠ(2))) <: IntegratorVPRK

end
