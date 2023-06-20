using GeometricIntegrators.Integrators
using GeometricEquations.Tests.HarmonicOscillator

using GeometricIntegrators.Methods: VPRK, VPRKMethod

iode = iodeproblem()
lode = lodeproblem()


@testset "$(rpad("Variational Partitioned Runge-Kutta Methods",80))" begin

    # @test VPRK(SRK3()) == VPRK(VPSRK3())

    # @test VPRK(Gauss(1)) == VPRK(VPRKGauss(1))
    @test VPRK(RadauIIA(2)) == VPRK(VPRKRadauIIA(2))
    @test VPRK(RadauIIB(2)) == VPRK(VPRKRadauIIB(2))

    @test VPRK(VPRKLobattoIII(2))  == VPRK(VPRKLobattoIII(2))
    @test VPRK(VPRKLobattoIIIA(2)) == VPRK(VPRKLobattoIIIA(2))
    @test VPRK(VPRKLobattoIIIB(2)) == VPRK(VPRKLobattoIIIB(2))
    @test VPRK(VPRKLobattoIIIC(2)) == VPRK(VPRKLobattoIIIC(2))
    @test VPRK(VPRKLobattoIIID(2)) == VPRK(VPRKLobattoIIID(2))
    @test VPRK(VPRKLobattoIIIE(2)) == VPRK(VPRKLobattoIIIE(2))
    @test VPRK(VPRKLobattoIIIF(2)) == VPRK(VPRKLobattoIIIF(2))
    @test VPRK(VPRKLobattoIIIG(2)) == VPRK(VPRKLobattoIIIG(2))

    # @test VPRK(LobattoIIIAIIIB(2)) == VPRK(VPRKLobattoIIIAIIIB(2))
    # @test VPRK(LobattoIIIBIIIA(2)) == VPRK(VPRKLobattoIIIBIIIA(2))
    # @test VPRK(LobattoIIIAIIIĀ(2)) == VPRK(VPRKLobattoIIIAIIIĀ(2))
    # @test VPRK(LobattoIIIBIIIB̄(2)) == VPRK(VPRKLobattoIIIBIIIB̄(2))
    @test VPRK(LobattoIIICIIIC̄(2)) == VPRK(VPRKLobattoIIICIIIC̄(2))
    @test VPRK(LobattoIIIDIIID̄(2)) == VPRK(VPRKLobattoIIIDIIID̄(2))
    @test VPRK(LobattoIIIEIIIĒ(2)) == VPRK(VPRKLobattoIIIEIIIĒ(2))
    @test VPRK(LobattoIIIFIIIF̄(2)) == VPRK(VPRKLobattoIIIFIIIF̄(2))
    @test VPRK(LobattoIIIGIIIḠ(2)) == VPRK(VPRKLobattoIIIGIIIḠ(2))


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

    @test typeof(VPRKLobattoIII(2)) <: VPRKMethod
    @test typeof(VPRKLobattoIII(3)) <: VPRKMethod
    @test typeof(VPRKLobattoIII(4)) <: VPRKMethod
    @test typeof(VPRKLobattoIII(5)) <: VPRKMethod

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


    @test typeof(Integrator(iode, VPSRK3())) <: IntegratorVPRK
    @test typeof(Integrator(iode, VPRKGauss(2))) <: IntegratorVPRK

    @test typeof(Integrator(iode, VPRKLobattoIII(2))) <: IntegratorVPRK
    @test typeof(Integrator(iode, VPRKLobattoIIIA(2))) <: IntegratorVPRK
    @test typeof(Integrator(iode, VPRKLobattoIIIB(2))) <: IntegratorVPRK
    @test typeof(Integrator(iode, VPRKLobattoIIIC(2))) <: IntegratorVPRK
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


@testset "$(rpad("Variational Partitioned Runge-Kutta Methods for Degenerate Lagrangians",80))" begin

    @test typeof(DegenerateVPRK(Gauss(1))) <: VPRKMethod
    @test typeof(DegenerateVPRK(Gauss(2))) <: VPRKMethod

    @test typeof(DegenerateVPRK(VPRKGauss(1))) <: VPRKMethod
    @test typeof(DegenerateVPRK(VPRKGauss(2))) <: VPRKMethod

end
