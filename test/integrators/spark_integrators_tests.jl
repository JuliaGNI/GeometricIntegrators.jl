
using GeometricIntegrators.Config
using GeometricIntegrators.Integrators
using GeometricIntegrators.Integrators.SPARK
using GeometricIntegrators.Tableaus
using GeometricIntegrators.Utils
using GeometricProblems.LotkaVolterra2d
using SimpleSolvers
using Test

SimpleSolvers.set_config(:nls_atol, 8eps())
SimpleSolvers.set_config(:nls_rtol, 2eps())

SimpleSolvers.set_config(:nls_atol_break, Inf)
SimpleSolvers.set_config(:nls_rtol_break, Inf)
SimpleSolvers.set_config(:nls_stol_break, Inf)

const Δt = 0.01
const nt = 10
const q₀ = [1.0, 1.0]
const parameters = (a₁=1.0, a₂=1.0, b₁=-1.0, b₂=-2.0)

ode  = lotka_volterra_2d_ode(q₀; params=parameters)
idae = lotka_volterra_2d_idae(q₀; params=parameters)
pdae = lotka_volterra_2d_pdae(q₀; params=parameters)
vdae = lotka_volterra_2d_vdae(q₀; params=parameters)
vdae_slrk = lotka_volterra_2d_slrk(q₀; params=parameters)

int  = IntegratorFIRK(ode, TableauGLRK(8), Δt)
sol  = integrate(ode, int, nt)
refx = sol.q[:,end]



@testset "$(rpad("SLRK integrators",80))" begin

    int = IntegratorSLRK(vdae_slrk, TableauSLRKLobattoIIIAB(2), Δt)
    sol = integrate(vdae_slrk, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = IntegratorSLRK(vdae_slrk, TableauSLRKLobattoIIIAB(3), Δt)
    sol = integrate(vdae_slrk, int, nt)
    @test rel_err(sol.q, refx) < 2E-11

    int = IntegratorSLRK(vdae_slrk, TableauSLRKLobattoIIIAB(4), Δt)
    sol = integrate(vdae_slrk, int, nt)
    @test rel_err(sol.q, refx) < 1E-15


    int = IntegratorSLRK(vdae_slrk, TableauSLRKLobattoIIIBA(2), Δt)
    sol = integrate(vdae_slrk, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = IntegratorSLRK(vdae_slrk, TableauSLRKLobattoIIIBA(3), Δt)
    sol = integrate(vdae_slrk, int, nt)
    @test rel_err(sol.q, refx) < 2E-11

    int = IntegratorSLRK(vdae_slrk, TableauSLRKLobattoIIIBA(4), Δt)
    sol = integrate(vdae_slrk, int, nt)
    @test rel_err(sol.q, refx) < 1E-15


    int = IntegratorSLRK(vdae_slrk, TableauSLRKLobattoIIICC̄(2), Δt)
    sol = integrate(vdae_slrk, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = IntegratorSLRK(vdae_slrk, TableauSLRKLobattoIIICC̄(3), Δt)
    sol = integrate(vdae_slrk, int, nt)
    @test rel_err(sol.q, refx) < 2E-11

    int = IntegratorSLRK(vdae_slrk, TableauSLRKLobattoIIICC̄(4), Δt)
    sol = integrate(vdae_slrk, int, nt)
    @test rel_err(sol.q, refx) < 1E-15


    int = IntegratorSLRK(vdae_slrk, TableauSLRKLobattoIIIC̄C(2), Δt)
    sol = integrate(vdae_slrk, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = IntegratorSLRK(vdae_slrk, TableauSLRKLobattoIIIC̄C(3), Δt)
    sol = integrate(vdae_slrk, int, nt)
    @test rel_err(sol.q, refx) < 2E-11

    int = IntegratorSLRK(vdae_slrk, TableauSLRKLobattoIIIC̄C(4), Δt)
    sol = integrate(vdae_slrk, int, nt)
    @test rel_err(sol.q, refx) < 1E-15


    int = IntegratorSLRK(vdae_slrk, TableauSLRKLobattoIIID(2), Δt)
    sol = integrate(vdae_slrk, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = IntegratorSLRK(vdae_slrk, TableauSLRKLobattoIIID(3), Δt)
    sol = integrate(vdae_slrk, int, nt)
    @test rel_err(sol.q, refx) < 2E-11

    int = IntegratorSLRK(vdae_slrk, TableauSLRKLobattoIIID(4), Δt)
    sol = integrate(vdae_slrk, int, nt)
    @test rel_err(sol.q, refx) < 2E-15


    int = IntegratorSLRK(vdae_slrk, TableauSLRKLobattoIIIE(2), Δt)
    sol = integrate(vdae_slrk, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = IntegratorSLRK(vdae_slrk, TableauSLRKLobattoIIIE(3), Δt)
    sol = integrate(vdae_slrk, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = IntegratorSLRK(vdae_slrk, TableauSLRKLobattoIIIE(4), Δt)
    sol = integrate(vdae_slrk, int, nt)
    @test rel_err(sol.q, refx) < 1E-15

end


@testset "$(rpad("SPARK integrators",80))" begin

    int = Integrator(idae, TableauSPARKGLRK(1), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = Integrator(idae, TableauSPARKGLRK(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = Integrator(idae, TableauSPARKLobattoIIIAIIIB(3), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = Integrator(idae, TableauSPARKLobattoIIIAIIIB(4), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 2E-10

    int = Integrator(idae, TableauSPARKGLRKLobattoIIIAIIIB(1), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 5E-4

    int = Integrator(idae, TableauSPARKGLRKLobattoIIIAIIIB(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 3E-4

    int = Integrator(idae, TableauSPARKGLRKLobattoIIIAIIIB(3), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 2E-4

    int = Integrator(idae, TableauSPARKGLRKLobattoIIIBIIIA(1), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 5E-4

    int = Integrator(idae, TableauSPARKGLRKLobattoIIIBIIIA(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 3E-4

    int = Integrator(idae, TableauSPARKGLRKLobattoIIIBIIIA(3), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 2E-4

end


@testset "$(rpad("VPARK integrators",80))" begin

    int = Integrator(idae, TableauSymplecticProjection(:glrk1ps, CoefficientsGLRK(1), CoefficientsGLRK(1)), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = Integrator(idae, TableauSymplecticProjection(:glrk2ps, CoefficientsGLRK(2), CoefficientsGLRK(2)), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = Integrator(idae, TableauGLRKpSymplectic(1), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = Integrator(idae, TableauGLRKpSymplectic(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = Integrator(idae, TableauLobattoIIIAIIIBpSymplectic(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 2E-6

    int = Integrator(idae, TableauLobattoIIIAIIIBpSymplectic(3), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 2E-11

    int = Integrator(idae, TableauLobattoIIIAIIIBpSymplectic(4), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 4E-16

    int = Integrator(idae, TableauLobattoIIIBIIIApSymplectic(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 2E-6

    int = Integrator(idae, TableauLobattoIIIBIIIApSymplectic(3), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 4E-11

    int = Integrator(idae, TableauLobattoIIIBIIIApSymplectic(4), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 4E-16

end


@testset "$(rpad("VSPARK integrators",80))" begin

    ### VSPARK Integrators ###

    int = IntegratorVSPARK(idae, TableauSPARKGLRKLobattoIIIAIIIB(1), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = IntegratorVSPARK(idae, TableauSPARKGLRKLobattoIIIAIIIB(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = IntegratorVSPARK(idae, TableauSPARKGLRKLobattoIIIAIIIB(3), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 5E-16

    int = IntegratorVSPARK(idae, TableauSPARKGLRKLobattoIIIBIIIA(1), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = IntegratorVSPARK(idae, TableauSPARKGLRKLobattoIIIBIIIA(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = IntegratorVSPARK(idae, TableauSPARKGLRKLobattoIIIBIIIA(3), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 5E-16


    ### VSPARKprimary Integrators ###

    int = Integrator(idae, TableauVSPARKGLRKpMidpoint(1), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = Integrator(idae, TableauVSPARKGLRKpMidpoint(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = Integrator(idae, TableauVSPARKGLRKpSymplectic(1), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = Integrator(idae, TableauVSPARKGLRKpSymplectic(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = Integrator(idae, TableauVSPARKGLRKpSymmetric(1), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = Integrator(idae, TableauVSPARKGLRKpSymmetric(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = IntegratorVSPARKprimary(idae, TableauVSPARKLobattoIIIAIIIBpSymmetric(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 2E-6

    int = IntegratorVSPARKprimary(idae, TableauVSPARKLobattoIIIAIIIBpSymmetric(3), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 5E-11

    int = IntegratorVSPARKprimary(idae, TableauVSPARKLobattoIIIAIIIBpSymmetric(4), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 2E-15

    int = IntegratorVSPARKprimary(idae, TableauVSPARKLobattoIIIBIIIApSymmetric(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 2E-6

    # int = IntegratorVSPARKprimary(idae, TableauVSPARKLobattoIIIBIIIApSymmetric(3), Δt)
    # sol = integrate(idae, int, nt)
    # @test rel_err(sol.q, refx) < 5E-11

    int = IntegratorVSPARKprimary(idae, TableauVSPARKLobattoIIIBIIIApSymmetric(4), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 2E-15


    ### VSPARKsecondary Integrators ###

    int = Integrator(vdae, TableauVSPARKLobattoIIIAB(2), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(vdae, TableauVSPARKLobattoIIIAB(3), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 2E-11

    int = Integrator(vdae, TableauVSPARKLobattoIIIAB(4), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-15


    int = Integrator(vdae, TableauVSPARKLobattoIIIBA(2), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(vdae, TableauVSPARKLobattoIIIBA(3), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 2E-11

    int = Integrator(vdae, TableauVSPARKLobattoIIIBA(4), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-15


    int = Integrator(vdae, TableauVSPARKLobattoIIICC̄(2), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(vdae, TableauVSPARKLobattoIIICC̄(3), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 2E-11

    int = Integrator(vdae, TableauVSPARKLobattoIIICC̄(4), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-15


    int = Integrator(vdae, TableauVSPARKLobattoIIIC̄C(2), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(vdae, TableauVSPARKLobattoIIIC̄C(3), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 2E-11

    int = Integrator(vdae, TableauVSPARKLobattoIIIC̄C(4), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-15


    int = Integrator(vdae, TableauVSPARKLobattoIIID(2), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(vdae, TableauVSPARKLobattoIIID(3), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 2E-11

    int = Integrator(vdae, TableauVSPARKLobattoIIID(4), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 2E-15


    int = Integrator(vdae, TableauVSPARKLobattoIIIE(2), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = Integrator(vdae, TableauVSPARKLobattoIIIE(3), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = Integrator(vdae, TableauVSPARKLobattoIIIE(4), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-15


    int = Integrator(vdae, TableauVSPARKGLRKLobattoIIIAB(1), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(vdae, TableauVSPARKGLRKLobattoIIIAB(2), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = Integrator(vdae, TableauVSPARKGLRKLobattoIIIAB(3), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-15


    int = Integrator(vdae, TableauVSPARKGLRKLobattoIIIBA(1), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(vdae, TableauVSPARKGLRKLobattoIIIBA(2), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = Integrator(vdae, TableauVSPARKGLRKLobattoIIIBA(3), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-15


    int = Integrator(vdae, TableauVSPARKGLRKLobattoIIICC̄(1), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(vdae, TableauVSPARKGLRKLobattoIIICC̄(2), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = Integrator(vdae, TableauVSPARKGLRKLobattoIIICC̄(3), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 2E-15


    int = Integrator(vdae, TableauVSPARKGLRKLobattoIIIC̄C(1), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(vdae, TableauVSPARKGLRKLobattoIIIC̄C(2), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = Integrator(vdae, TableauVSPARKGLRKLobattoIIIC̄C(3), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 2E-15


    int = Integrator(vdae, TableauVSPARKGLRKLobattoIIID(1), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(vdae, TableauVSPARKGLRKLobattoIIID(2), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = Integrator(vdae, TableauVSPARKGLRKLobattoIIID(3), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-15


    int = Integrator(vdae, TableauVSPARKGLRKLobattoIIIE(1), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = Integrator(vdae, TableauVSPARKGLRKLobattoIIIE(2), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = Integrator(vdae, TableauVSPARKGLRKLobattoIIIE(3), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-15
end


@testset "$(rpad("HPARK integrators",80))" begin

    int = Integrator(pdae, TableauHPARKGLRK(1), Δt)
    sol = integrate(pdae, int, nt)
    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 2E-6

    int = Integrator(pdae, TableauHPARKGLRK(2), Δt)
    sol = integrate(pdae, int, nt)
    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 8E-7

    # int = Integrator(pdae, TableauHPARKLobattoIIIAIIIB(2), Δt)
    # sol = integrate(pdae, int, nt)
    # TODO
    # println(rel_err(sol.q, refx))
    # @test rel_err(sol.q, refx) < 2E-2

    int = Integrator(pdae, TableauHPARKLobattoIIIAIIIB(3), Δt)
    sol = integrate(pdae, int, nt)
    # TODO
    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 8E-2

    int = Integrator(pdae, TableauHPARKLobattoIIIAIIIB(4), Δt)
    sol = integrate(pdae, int, nt)
    # TODO
    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 8E-4

    # int = Integrator(pdae, TableauHPARKLobattoIIIBIIIA(2), Δt)
    # sol = integrate(pdae, int, nt)
    # TODO
    # println(rel_err(sol.q, refx))
    # @test rel_err(sol.q, refx) < 2E-2

    int = Integrator(pdae, TableauHPARKLobattoIIIBIIIA(3), Δt)
    sol = integrate(pdae, int, nt)
    # TODO
    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 4E-3

    int = Integrator(pdae, TableauHPARKLobattoIIIBIIIA(4), Δt)
    sol = integrate(pdae, int, nt)
    # TODO
    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 2E-2

end


@testset "$(rpad("HSPARK integrators",80))" begin
    ### HSPARK Integrators ###

    int = IntegratorHSPARK(pdae, TableauSPARKGLRK(1), Δt)
    sol = integrate(pdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = IntegratorHSPARK(pdae, TableauSPARKGLRK(2), Δt)
    sol = integrate(pdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11


    ### HSPARKprimary Integrators ###

    int = Integrator(pdae, TableauHSPARKGLRKpSymmetric(1), Δt)
    sol = integrate(pdae, int, nt)
    # TODO
    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(pdae, TableauHSPARKGLRKpSymmetric(2), Δt)
    sol = integrate(pdae, int, nt)
    # TODO
    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(pdae, TableauHSPARKLobattoIIIAIIIBpSymmetric(2), Δt)
    sol = integrate(pdae, int, nt)
    # TODO
    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(pdae, TableauHSPARKLobattoIIIAIIIBpSymmetric(3), Δt)
    sol = integrate(pdae, int, nt)
    # TODO
    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(pdae, TableauHSPARKLobattoIIIAIIIBpSymmetric(4), Δt)
    sol = integrate(pdae, int, nt)
    # TODO
    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(pdae, TableauHSPARKLobattoIIIBIIIApSymmetric(2), Δt)
    sol = integrate(pdae, int, nt)
    # TODO
    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(pdae, TableauHSPARKLobattoIIIBIIIApSymmetric(3), Δt)
    sol = integrate(pdae, int, nt)
    # TODO
    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(pdae, TableauHSPARKLobattoIIIBIIIApSymmetric(4), Δt)
    sol = integrate(pdae, int, nt)
    # TODO
    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 4E-6


    ### HSPARKsecondary Integrators ###

    # int = IntegratorHSPARKsecondary(vdae, TableauHSPARKLobattoIIIAB(2), Δt)
    # sol = integrate(vdae, int, nt)
    # println(rel_err(sol.q, refx))
    # @test rel_err(sol.q, refx) < 4E-6

    # int = IntegratorHSPARKsecondary(vdae, TableauHSPARKLobattoIIIAB(3), Δt)
    # sol = integrate(vdae, int, nt)
    # println(rel_err(sol.q, refx))
    # @test rel_err(sol.q, refx) < 2E-11

    # int = IntegratorHSPARKsecondary(vdae, TableauHSPARKLobattoIIIAB(4), Δt)
    # sol = integrate(vdae, int, nt)
    # println(rel_err(sol.q, refx))
    # @test rel_err(sol.q, refx) < 1E-15

end
