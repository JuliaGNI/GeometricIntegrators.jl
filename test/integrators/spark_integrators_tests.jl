
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


@testset "$(rpad("VPARK integrators",80))" begin

    int = Integrator(idae, getTableauSymplecticProjection(:pglrk1ps, getCoefficientsGLRK(1), getCoefficientsGLRK(1)), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = Integrator(idae, getTableauSymplecticProjection(:pglrk2ps, getCoefficientsGLRK(2), getCoefficientsGLRK(2)), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = Integrator(idae, getTableauGLRKpSymplectic(1), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = Integrator(idae, getTableauGLRKpSymplectic(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = Integrator(idae, getTableauLobIIIAIIIB2pSymplectic(), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 2E-6

    int = Integrator(idae, getTableauLobIIIAIIIB3pSymplectic(), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 2E-11

    int = Integrator(idae, getTableauLobIIIAIIIB4pSymplectic(), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 4E-16

end


@testset "$(rpad("VSPARK integrators",80))" begin

    ### VSPARK Integrators ###

    int = IntegratorVSPARK(idae, getTableauSPARKGLRKLobIIIAIIIB(1), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = IntegratorVSPARK(idae, getTableauSPARKGLRKLobIIIAIIIB(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = IntegratorVSPARK(idae, getTableauSPARKGLRKLobIIIAIIIB(3), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 4E-16


    ### VSPARKprimary Integrators ###

    int = Integrator(idae, getTableauVSPARKGLRKpMidpoint(1), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = Integrator(idae, getTableauVSPARKGLRKpMidpoint(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = Integrator(idae, getTableauVSPARKGLRKpSymplectic(1), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = Integrator(idae, getTableauVSPARKGLRKpSymplectic(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = Integrator(idae, getTableauVSPARKGLRKpSymmetric(1), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = Integrator(idae, getTableauVSPARKGLRKpSymmetric(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = IntegratorVSPARKprimary(idae, getTableauVSPARKLobIIIAIIIB2pSymmetric(), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 2E-6

    int = IntegratorVSPARKprimary(idae, getTableauVSPARKLobIIIAIIIB3pSymmetric(), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 2E-11

    int = IntegratorVSPARKprimary(idae, getTableauVSPARKLobIIIAIIIB4pSymmetric(), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, refx) < 8E-16


    ### VSPARKsecondary Integrators ###

    int = Integrator(vdae, getTableauVSPARKLobIIIAB(2), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(vdae, getTableauVSPARKLobIIIAB(3), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 2E-11

    int = Integrator(vdae, getTableauVSPARKLobIIIAB(4), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-15


    int = Integrator(vdae, getTableauVSPARKLobIIIC(2), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(vdae, getTableauVSPARKLobIIIC(3), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 2E-11

    int = Integrator(vdae, getTableauVSPARKLobIIIC(4), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-15


    int = Integrator(vdae, getTableauVSPARKLobIIID(2), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(vdae, getTableauVSPARKLobIIID(3), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 2E-11

    int = Integrator(vdae, getTableauVSPARKLobIIID(4), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 2E-15


    int = Integrator(vdae, getTableauVSPARKLobIIIE(2), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = Integrator(vdae, getTableauVSPARKLobIIIE(3), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = Integrator(vdae, getTableauVSPARKLobIIIE(4), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-15


    int = Integrator(vdae, getTableauVSPARKGLRKLobIIIAB(1), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(vdae, getTableauVSPARKGLRKLobIIIAB(2), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = Integrator(vdae, getTableauVSPARKGLRKLobIIIAB(3), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-15


    int = Integrator(vdae, getTableauVSPARKGLRKLobIIIC(1), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(vdae, getTableauVSPARKGLRKLobIIIC(2), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = Integrator(vdae, getTableauVSPARKGLRKLobIIIC(3), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 2E-15


    int = Integrator(vdae, getTableauVSPARKGLRKLobIIID(1), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(vdae, getTableauVSPARKGLRKLobIIID(2), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = Integrator(vdae, getTableauVSPARKGLRKLobIIID(3), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-15


    int = Integrator(vdae, getTableauVSPARKGLRKLobIIIE(1), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = Integrator(vdae, getTableauVSPARKGLRKLobIIIE(2), Δt)
    sol = integrate(vdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11

    int = Integrator(vdae, getTableauVSPARKGLRKLobIIIE(3), Δt)
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

    int = Integrator(pdae, getTableauHPARKGLRK(1), Δt)
    sol = integrate(pdae, int, nt)

    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 2E-6

    int = Integrator(pdae, getTableauHPARKGLRK(2), Δt)
    sol = integrate(pdae, int, nt)

    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 8E-7

    int = Integrator(pdae, getTableauHPARKLobIIIAIIIB2(), Δt)
    # sol = integrate(pdae, int, nt)
    # TODO
    # println(rel_err(sol.q, refx))
    # @test rel_err(sol.q, refx) < 2E-2

    int = Integrator(pdae, getTableauHPARKLobIIIAIIIB3(), Δt)
    # sol = integrate(pdae, int, nt)
    # TODO
    # println(rel_err(sol.q, refx))
    # @test rel_err(sol.q, refx) < 8E-2

    int = Integrator(pdae, getTableauHPARKLobIIIAIIIB4(), Δt)
    # sol = integrate(pdae, int, nt)
    # TODO
    # println(rel_err(sol.q, refx))
    # @test rel_err(sol.q, refx) < 8E-2
end


@testset "$(rpad("HSPARK integrators",80))" begin
    ### HSPARK Integrators ###

    int = IntegratorHSPARK(pdae, getTableauSPARKGLRK(1), Δt)
    sol = integrate(pdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-6

    int = IntegratorHSPARK(pdae, getTableauSPARKGLRK(2), Δt)
    sol = integrate(pdae, int, nt)
    @test rel_err(sol.q, refx) < 1E-11


    ### HSPARKprimary Integrators ###

    int = Integrator(pdae, getTableauHSPARKGLRKpSymmetric(1), Δt)
    sol = integrate(pdae, int, nt)

    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(pdae, getTableauHSPARKGLRKpSymmetric(2), Δt)
    sol = integrate(pdae, int, nt)

    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(pdae, getTableauHSPARKLobIIIAIIIB2pSymmetric(), Δt)
    sol = integrate(pdae, int, nt)

    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(pdae, getTableauHSPARKLobIIIAIIIB3pSymmetric(), Δt)
    sol = integrate(pdae, int, nt)

    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 4E-6

    int = Integrator(pdae, getTableauHSPARKLobIIIAIIIB4pSymmetric(), Δt)
    sol = integrate(pdae, int, nt)

    # println(rel_err(sol.q, refx))
    @test rel_err(sol.q, refx) < 4E-6


    ### HSPARKsecondary Integrators ###

    # int = IntegratorHSPARKsecondary(vdae, getTableauHSPARKLobIIIAB(2), Δt)
    # sol = integrate(vdae, int, nt)
    # println(rel_err(sol.q, refx))
    # @test rel_err(sol.q, refx) < 4E-6

    # int = IntegratorHSPARKsecondary(vdae, getTableauHSPARKLobIIIAB(3), Δt)
    # sol = integrate(vdae, int, nt)
    # println(rel_err(sol.q, refx))
    # @test rel_err(sol.q, refx) < 2E-11

    # int = IntegratorHSPARKsecondary(vdae, getTableauHSPARKLobIIIAB(4), Δt)
    # sol = integrate(vdae, int, nt)
    # println(rel_err(sol.q, refx))
    # @test rel_err(sol.q, refx) < 1E-15

end
