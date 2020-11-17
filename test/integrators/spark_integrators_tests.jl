
using GeometricIntegrators.Config
using GeometricIntegrators.Integrators
using GeometricIntegrators.Integrators.SPARK
using GeometricIntegrators.Solvers
using GeometricIntegrators.Tableaus
using GeometricIntegrators.TestProblems.LotkaVolterra2dProblem
using GeometricIntegrators.Utils
using Test

using GeometricIntegrators.TestProblems.LotkaVolterra2dProblem: Δt, nt

set_config(:nls_atol, 8eps())
set_config(:nls_rtol, 2eps())

set_config(:nls_atol_break, Inf)
set_config(:nls_rtol_break, Inf)
set_config(:nls_stol_break, Inf)

ode  = lotka_volterra_2d_ode()
idae = lotka_volterra_2d_idae()
pdae = lotka_volterra_2d_pdae()
vdae = lotka_volterra_2d_vdae()
vdae_slrk = lotka_volterra_2d_vdae_slrk()

int  = IntegratorFIRK(ode, getTableauGLRK(8), Δt)
sol  = integrate(ode, int, nt)
refx = sol.q[:,end]



@testset "$(rpad("SLRK integrators",80))" begin

    int = IntegratorSLRK(vdae_slrk, getTableauSLRKLobIIIAB(2), Δt)
    sol = integrate(vdae_slrk, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = IntegratorSLRK(vdae_slrk, getTableauSLRKLobIIIAB(3), Δt)
    sol = integrate(vdae_slrk, int, nt)
    @test rel_err(sol.q, refx) < 2E-11

    int = IntegratorSLRK(vdae_slrk, getTableauSLRKLobIIIAB(4), Δt)
    sol = integrate(vdae_slrk, int, nt)
    @test rel_err(sol.q, refx) < 1E-15


    int = IntegratorSLRK(vdae_slrk, getTableauSLRKLobIIIC(2), Δt)
    sol = integrate(vdae_slrk, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = IntegratorSLRK(vdae_slrk, getTableauSLRKLobIIIC(3), Δt)
    sol = integrate(vdae_slrk, int, nt)
    @test rel_err(sol.q, refx) < 2E-11

    int = IntegratorSLRK(vdae_slrk, getTableauSLRKLobIIIC(4), Δt)
    sol = integrate(vdae_slrk, int, nt)
    @test rel_err(sol.q, refx) < 1E-15


    int = IntegratorSLRK(vdae_slrk, getTableauSLRKLobIIID(2), Δt)
    sol = integrate(vdae_slrk, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = IntegratorSLRK(vdae_slrk, getTableauSLRKLobIIID(3), Δt)
    sol = integrate(vdae_slrk, int, nt)
    @test rel_err(sol.q, refx) < 2E-11

    int = IntegratorSLRK(vdae_slrk, getTableauSLRKLobIIID(4), Δt)
    sol = integrate(vdae_slrk, int, nt)
    @test rel_err(sol.q, refx) < 2E-15


    int = IntegratorSLRK(vdae_slrk, getTableauSLRKLobIIIE(2), Δt)
    sol = integrate(vdae_slrk, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = IntegratorSLRK(vdae_slrk, getTableauSLRKLobIIIE(3), Δt)
    sol = integrate(vdae_slrk, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = IntegratorSLRK(vdae_slrk, getTableauSLRKLobIIIE(4), Δt)
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

    int = Integrator(idae, TableauSPARKLobIIIAIIIB(3), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = Integrator(idae, TableauSPARKLobIIIAIIIB(4), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 2E-10

    int = Integrator(idae, TableauSPARKGLRKLobIIIAIIIB(1), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 5E-4

    int = Integrator(idae, TableauSPARKGLRKLobIIIAIIIB(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 3E-4

    int = Integrator(idae, TableauSPARKGLRKLobIIIAIIIB(3), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 2E-4

    int = Integrator(idae, TableauSPARKGLRKLobIIIBIIIA(1), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 5E-4

    int = Integrator(idae, TableauSPARKGLRKLobIIIBIIIA(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 3E-4

    int = Integrator(idae, TableauSPARKGLRKLobIIIBIIIA(3), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 2E-4

end


@testset "$(rpad("VPARK integrators",80))" begin

    int = Integrator(idae, TableauSymplecticProjection(:glrk1ps, getCoefficientsGLRK(1), getCoefficientsGLRK(1)), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = Integrator(idae, TableauSymplecticProjection(:glrk2ps, getCoefficientsGLRK(2), getCoefficientsGLRK(2)), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = Integrator(idae, TableauGLRKpSymplectic(1), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = Integrator(idae, TableauGLRKpSymplectic(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = Integrator(idae, TableauLobIIIAIIIBpSymplectic(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 2E-6

    int = Integrator(idae, TableauLobIIIAIIIBpSymplectic(3), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 2E-11

    int = Integrator(idae, TableauLobIIIAIIIBpSymplectic(4), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 4E-16

    int = Integrator(idae, TableauLobIIIBIIIApSymplectic(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 2E-6

    int = Integrator(idae, TableauLobIIIBIIIApSymplectic(3), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 4E-11

    int = Integrator(idae, TableauLobIIIBIIIApSymplectic(4), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 4E-16

end


@testset "$(rpad("VSPARK integrators",80))" begin

    ### VSPARK Integrators ###

    int = IntegratorVSPARK(idae, TableauSPARKGLRKLobIIIAIIIB(1), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = IntegratorVSPARK(idae, TableauSPARKGLRKLobIIIAIIIB(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = IntegratorVSPARK(idae, TableauSPARKGLRKLobIIIAIIIB(3), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 5E-16

    int = IntegratorVSPARK(idae, TableauSPARKGLRKLobIIIBIIIA(1), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = IntegratorVSPARK(idae, TableauSPARKGLRKLobIIIBIIIA(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = IntegratorVSPARK(idae, TableauSPARKGLRKLobIIIBIIIA(3), Δt)
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

    int = IntegratorVSPARKprimary(idae, TableauVSPARKLobIIIAIIIBpSymmetric(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 2E-6

    int = IntegratorVSPARKprimary(idae, TableauVSPARKLobIIIAIIIBpSymmetric(3), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 5E-11

    int = IntegratorVSPARKprimary(idae, TableauVSPARKLobIIIAIIIBpSymmetric(4), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 2E-15

    int = IntegratorVSPARKprimary(idae, TableauVSPARKLobIIIBIIIApSymmetric(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 2E-6

    # int = IntegratorVSPARKprimary(idae, TableauVSPARKLobIIIBIIIApSymmetric(3), Δt)
    # sol = integrate(idae, int, nt)
    # @test rel_err(sol.q, refx) < 5E-11

    int = IntegratorVSPARKprimary(idae, TableauVSPARKLobIIIBIIIApSymmetric(4), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 2E-15


    ### VSPARKsecondary Integrators ###

    int = Integrator(vdae, TableauVSPARKLobIIIAB(2), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(vdae, TableauVSPARKLobIIIAB(3), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 2E-11

    int = Integrator(vdae, TableauVSPARKLobIIIAB(4), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-15


    int = Integrator(vdae, TableauVSPARKLobIIIBA(2), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(vdae, TableauVSPARKLobIIIBA(3), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 2E-11

    int = Integrator(vdae, TableauVSPARKLobIIIBA(4), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-15


    int = Integrator(vdae, TableauVSPARKLobIIIC(2), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(vdae, TableauVSPARKLobIIIC(3), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 2E-11

    int = Integrator(vdae, TableauVSPARKLobIIIC(4), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-15


    int = Integrator(vdae, TableauVSPARKLobIIID(2), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(vdae, TableauVSPARKLobIIID(3), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 2E-11

    int = Integrator(vdae, TableauVSPARKLobIIID(4), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 2E-15


    int = Integrator(vdae, TableauVSPARKLobIIIE(2), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = Integrator(vdae, TableauVSPARKLobIIIE(3), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = Integrator(vdae, TableauVSPARKLobIIIE(4), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-15


    int = Integrator(vdae, TableauVSPARKGLRKLobIIIAB(1), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(vdae, TableauVSPARKGLRKLobIIIAB(2), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = Integrator(vdae, TableauVSPARKGLRKLobIIIAB(3), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-15


    int = Integrator(vdae, TableauVSPARKGLRKLobIIIBA(1), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(vdae, TableauVSPARKGLRKLobIIIBA(2), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = Integrator(vdae, TableauVSPARKGLRKLobIIIBA(3), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-15


    int = Integrator(vdae, TableauVSPARKGLRKLobIIIC(1), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(vdae, TableauVSPARKGLRKLobIIIC(2), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = Integrator(vdae, TableauVSPARKGLRKLobIIIC(3), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 2E-15


    int = Integrator(vdae, TableauVSPARKGLRKLobIIID(1), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(vdae, TableauVSPARKGLRKLobIIID(2), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = Integrator(vdae, TableauVSPARKGLRKLobIIID(3), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-15


    int = Integrator(vdae, TableauVSPARKGLRKLobIIIE(1), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = Integrator(vdae, TableauVSPARKGLRKLobIIIE(2), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = Integrator(vdae, TableauVSPARKGLRKLobIIIE(3), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-15
end


@testset "$(rpad("HPARK integrators",80))" begin

    int = Integrator(pdae, getTableauHPARK(:hpark_glrk1, getCoefficientsGLRK(1), getCoefficientsGLRK(1)), Δt)
    sol = integrate(pdae, int, nt)
    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 2E-6

    int = Integrator(pdae, getTableauHPARK(:hpark_glrk2, getCoefficientsGLRK(2), getCoefficientsGLRK(2)), Δt)
    sol = integrate(pdae, int, nt)
    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 8E-7

    int = Integrator(pdae, TableauHPARKGLRK(1), Δt)
    sol = integrate(pdae, int, nt)
    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 2E-6

    int = Integrator(pdae, TableauHPARKGLRK(2), Δt)
    sol = integrate(pdae, int, nt)
    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 8E-7

    # int = Integrator(pdae, TableauHPARKLobIIIAIIIB(2), Δt)
    # sol = integrate(pdae, int, nt)
    # TODO
    # println(rel_err(sol.q, refx))
    # @test rel_err(sol.q, refx) < 2E-2

    int = Integrator(pdae, TableauHPARKLobIIIAIIIB(3), Δt)
    sol = integrate(pdae, int, nt)
    # TODO
    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 8E-2

    int = Integrator(pdae, TableauHPARKLobIIIAIIIB(4), Δt)
    sol = integrate(pdae, int, nt)
    # TODO
    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 8E-4

    # int = Integrator(pdae, TableauHPARKLobIIIBIIIA(2), Δt)
    # sol = integrate(pdae, int, nt)
    # TODO
    # println(rel_err(sol.q, refx))
    # @test rel_err(sol.q, refx) < 2E-2

    int = Integrator(pdae, TableauHPARKLobIIIBIIIA(3), Δt)
    sol = integrate(pdae, int, nt)
    # TODO
    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 4E-3

    int = Integrator(pdae, TableauHPARKLobIIIBIIIA(4), Δt)
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

    int = Integrator(pdae, TableauHSPARKLobIIIAIIIBpSymmetric(2), Δt)
    sol = integrate(pdae, int, nt)
    # TODO
    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(pdae, TableauHSPARKLobIIIAIIIBpSymmetric(3), Δt)
    sol = integrate(pdae, int, nt)
    # TODO
    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(pdae, TableauHSPARKLobIIIAIIIBpSymmetric(4), Δt)
    sol = integrate(pdae, int, nt)
    # TODO
    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(pdae, TableauHSPARKLobIIIBIIIApSymmetric(2), Δt)
    sol = integrate(pdae, int, nt)
    # TODO
    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(pdae, TableauHSPARKLobIIIBIIIApSymmetric(3), Δt)
    sol = integrate(pdae, int, nt)
    # TODO
    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(pdae, TableauHSPARKLobIIIBIIIApSymmetric(4), Δt)
    sol = integrate(pdae, int, nt)
    # TODO
    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 4E-6


    ### HSPARKsecondary Integrators ###

    # int = IntegratorHSPARKsecondary(vdae, TableauHSPARKLobIIIAB(2), Δt)
    # sol = integrate(vdae, int, nt)
    # println(rel_err(sol.q, refx))
    # @test rel_err(sol.q, refx) < 4E-6

    # int = IntegratorHSPARKsecondary(vdae, TableauHSPARKLobIIIAB(3), Δt)
    # sol = integrate(vdae, int, nt)
    # println(rel_err(sol.q, refx))
    # @test rel_err(sol.q, refx) < 2E-11

    # int = IntegratorHSPARKsecondary(vdae, TableauHSPARKLobIIIAB(4), Δt)
    # sol = integrate(vdae, int, nt)
    # println(rel_err(sol.q, refx))
    # @test rel_err(sol.q, refx) < 1E-15

end
