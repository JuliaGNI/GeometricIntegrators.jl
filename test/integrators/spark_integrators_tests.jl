
using GeometricIntegrators.Config
using GeometricIntegrators.Integrators
using GeometricIntegrators.Solvers
using GeometricIntegrators.Tableaus
using GeometricIntegrators.TestProblems.LotkaVolterra
using GeometricIntegrators.Utils
using Test

set_config(:nls_atol, 8eps())
set_config(:nls_rtol, 2eps())
set_config(:nls_solver, QuasiNewtonSolver)
set_config(:jacobian_autodiff, true)
# set_config(:nls_nmax, 3)

using GeometricIntegrators.TestProblems.LotkaVolterra: Δt, nt

# Δt = 0.001
# nt = 100

idae = lotka_volterra_2d_idae()
pdae = lotka_volterra_2d_pdae()
vdae = lotka_volterra_2d_vdae()

int  = IntegratorFIRK(lotka_volterra_2d_ode(), getTableauGLRK(8), Δt)
sol  = integrate(int, nt)
refx = sol.q[:,end]

# set_config(:nls_nmax, 1)

# @testset "$(rpad("VPARK integrators",80))" begin
#
#     dint = Integrator(idae, getTableauSymplecticProjection(:pglrk1ps, getCoefficientsGLRK(1), getCoefficientsGLRK(1)), Δt)
#     dsol = integrate(dint, nt)
#
#     # println(rel_err(dsol.q, refx))
#     @test rel_err(dsol.q, refx) < 1E-6
#
#     dint = Integrator(idae, getTableauSymplecticProjection(:pglrk2ps, getCoefficientsGLRK(2), getCoefficientsGLRK(2)), Δt)
#     dsol = integrate(dint, nt)
#
#     # println(rel_err(dsol.q, refx))
#     @test rel_err(dsol.q, refx) < 1E-11
#
#     dint = Integrator(idae, getTableauGLRKpSymplectic(1), Δt)
#     dsol = integrate(dint, nt)
#
#     # println(rel_err(dsol.q, refx))
#     @test rel_err(dsol.q, refx) < 1E-6
#
#     dint = Integrator(idae, getTableauGLRKpSymplectic(2), Δt)
#     dsol = integrate(dint, nt)
#
#     # println(rel_err(dsol.q, refx))
#     @test rel_err(dsol.q, refx) < 1E-11
#
#     dint = Integrator(idae, getTableauLobIIIAIIIB2pSymplectic(), Δt)
#     dsol = integrate(dint, nt)
#
#     # println(rel_err(dsol.q, refx))
#     @test rel_err(dsol.q, refx) < 2E-6
#
#     dint = Integrator(idae, getTableauLobIIIAIIIB3pSymplectic(), Δt)
#     dsol = integrate(dint, nt)
#
#     # println(rel_err(dsol.q, refx))
#     @test rel_err(dsol.q, refx) < 8E-5
#
# end


@testset "$(rpad("VSPARK integrators",80))" begin
    # ### VSPARK Integrators ###
    #
    # dint = IntegratorVSPARK(idae, getTableauSPARKGLRK(1), Δt)
    # # dsol = integrate(dint, nt)
    # # TODO
    # # println(rel_err(dsol.q, refx))
    # # @test rel_err(dsol.q, refx) < 1E-6
    #
    #
    # ### VSPARKprimary Integrators ###
    #
    # dint = Integrator(idae, getTableauVSPARKGLRKpMidpoint(1), Δt)
    # dsol = integrate(dint, nt)
    #
    # # println(rel_err(dsol.q, refx))
    # @test rel_err(dsol.q, refx) < 1E-6
    #
    # dint = Integrator(idae, getTableauVSPARKGLRKpMidpoint(2), Δt)
    # dsol = integrate(dint, nt)
    #
    # # println(rel_err(dsol.q, refx))
    # @test rel_err(dsol.q, refx) < 1E-11
    #
    # dint = Integrator(idae, getTableauVSPARKGLRKpSymplectic(1), Δt)
    # dsol = integrate(dint, nt)
    #
    # # println(rel_err(dsol.q, refx))
    # @test rel_err(dsol.q, refx) < 1E-6
    #
    # dint = Integrator(idae, getTableauVSPARKGLRKpSymplectic(2), Δt)
    # dsol = integrate(dint, nt)
    #
    # # println(rel_err(dsol.q, refx))
    # @test rel_err(dsol.q, refx) < 1E-11
    #
    # dint = Integrator(idae, getTableauVSPARKGLRKpSymmetric(1), Δt)
    # dsol = integrate(dint, nt)
    #
    # # println(rel_err(dsol.q, refx))
    # @test rel_err(dsol.q, refx) < 1E-6
    #
    # dint = Integrator(idae, getTableauVSPARKGLRKpSymmetric(2), Δt)
    # dsol = integrate(dint, nt)
    #
    # # println(rel_err(dsol.q, refx))
    # @test rel_err(dsol.q, refx) < 1E-11
    #
    # dint = IntegratorVSPARKprimary(idae, getTableauVSPARKLobIIIAIIIB2pSymmetric(), Δt)
    # dsol = integrate(dint, nt)
    #
    # # println(rel_err(dsol.q, refx))
    # @test rel_err(dsol.q, refx) < 2E-6
    #
    # dint = IntegratorVSPARKprimary(idae, getTableauVSPARKLobIIIAIIIB3pSymmetric(), Δt)
    # dsol = integrate(dint, nt)
    #
    # # println(rel_err(dsol.q, refx))
    # @test rel_err(dsol.q, refx) < 8E-5


    ### VSPARKsecondary Integrators ###

    # dint = Integrator(vdae, getTableauVSPARKGLRK(1), Δt)
    # dsol = integrate(dint, nt)
    #
    # println(rel_err(dsol.q, refx))
    # # @test rel_err(dsol.q, refx) < 4E-6
    #
    # dint = Integrator(vdae, getTableauVSPARKGLRK(2), Δt)
    # dsol = integrate(dint, nt)
    #
    # println(rel_err(dsol.q, refx))
    # # @test rel_err(dsol.q, refx) < 1E-11
    #
    # dint = Integrator(vdae, getTableauVSPARKGLRK(3), Δt)
    # dsol = integrate(dint, nt)
    #
    # println(rel_err(dsol.q, refx))
    # # @test rel_err(dsol.q, refx) < 1E-11

    # dint = Integrator(vdae, getTableauVSPARKLobIII(2), Δt)
    # dsol = integrate(dint, nt)
    #
    # println(rel_err(dsol.q, refx))
    # # @test rel_err(dsol.q, refx) < 4E-6
    #
    # dint = Integrator(vdae, getTableauVSPARKLobIII(3), Δt)
    # dsol = integrate(dint, nt)
    #
    # println(rel_err(dsol.q, refx))
    # # @test rel_err(dsol.q, refx) < 1E-11
    #
    # dint = Integrator(vdae, getTableauVSPARKLobIII(4), Δt)
    # dsol = integrate(dint, nt)
    #
    # println(rel_err(dsol.q, refx))
    # # @test rel_err(dsol.q, refx) < 1E-11

    println(getTableauVSPARKGLRKLobIII(1))
    dint = Integrator(vdae, getTableauVSPARKGLRKLobIII(1), Δt)
    dsol = integrate(dint, nt)

    println(rel_err(dsol.q, refx))
    # @test rel_err(dsol.q, refx) < 4E-6

    println(getTableauVSPARKGLRKLobIII(2))
    dint = Integrator(vdae, getTableauVSPARKGLRKLobIII(2), Δt)
    dsol = integrate(dint, nt)

    println(rel_err(dsol.q, refx))
    # @test rel_err(dsol.q, refx) < 1E-11
    #
    # dint = Integrator(vdae, getTableauVSPARKGLRKLobIII(3), Δt)
    # dsol = integrate(dint, nt)
    #
    # println(rel_err(dsol.q, refx))
    # # @test rel_err(dsol.q, refx) < 1E-11
    #
    # dint = Integrator(vdae, getTableauVSPARKGLRKLobIII(1,2), Δt)
    # dsol = integrate(dint, nt)
    #
    # println(rel_err(dsol.q, refx), "  (s=1,σ=2)")
    # # @test rel_err(dsol.q, refx) < 4E-6
    #
    # dint = Integrator(vdae, getTableauVSPARKGLRKLobIII(1,3), Δt)
    # dsol = integrate(dint, nt)
    #
    # println(rel_err(dsol.q, refx), "  (s=1,σ=3)")
    # # @test rel_err(dsol.q, refx) < 4E-6
    #
    # dint = Integrator(vdae, getTableauVSPARKGLRKLobIII(1,4), Δt)
    # dsol = integrate(dint, nt)
    #
    # println(rel_err(dsol.q, refx), "  (s=1,σ=4)")
    # # @test rel_err(dsol.q, refx) < 4E-6
    #
    # dint = Integrator(vdae, getTableauVSPARKGLRKLobIII(2,3), Δt)
    # dsol = integrate(dint, nt)
    #
    # println(rel_err(dsol.q, refx), "  (s=2,σ=3)")
    # # @test rel_err(dsol.q, refx) < 1E-11
    #
    # dint = Integrator(vdae, getTableauVSPARKGLRKLobIII(2,4), Δt)
    # dsol = integrate(dint, nt)
    #
    # println(rel_err(dsol.q, refx), "  (s=2,σ=4)")
    # # @test rel_err(dsol.q, refx) < 1E-11
    #
    # # dint = Integrator(vdae, getTableauVSPARKGLRKLobIII(3,2), Δt)
    # # dsol = integrate(dint, nt)
    # #
    # # println(rel_err(dsol.q, refx), "(s=3,σ=2)")
    # # # @test rel_err(dsol.q, refx) < 1E-11
    #
    # dint = Integrator(vdae, getTableauVSPARKGLRKLobIII(3,3), Δt)
    # dsol = integrate(dint, nt)
    #
    # println(rel_err(dsol.q, refx), "  (s=3,σ=3)")
    # # @test rel_err(dsol.q, refx) < 1E-11
    #
    # dint = Integrator(vdae, getTableauVSPARKGLRKLobIII(3,4), Δt)
    # dsol = integrate(dint, nt)
    #
    # println(rel_err(dsol.q, refx), "  (s=3,σ=4)")
    # # @test rel_err(dsol.q, refx) < 1E-11
    #
    # dint = Integrator(vdae, getTableauVSPARKGLRKLobIIIAB(1), Δt)
    # dsol = integrate(dint, nt)
    #
    # println(rel_err(dsol.q, refx))
    # # @test rel_err(dsol.q, refx) < 4E-6
    #
    # dint = Integrator(vdae, getTableauVSPARKGLRKLobIIIAB(2), Δt)
    # dsol = integrate(dint, nt)
    #
    # println(rel_err(dsol.q, refx))
    # # @test rel_err(dsol.q, refx) < 1E-11
    #
    # dint = Integrator(vdae, getTableauVSPARKGLRKLobIIIAB(3), Δt)
    # dsol = integrate(dint, nt)
    #
    # println(rel_err(dsol.q, refx))
    # # @test rel_err(dsol.q, refx) < 1E-11
    #
    # dint = Integrator(vdae, getTableauVSPARKGLRKLobIIIAC(1), Δt)
    # dsol = integrate(dint, nt)
    #
    # println(rel_err(dsol.q, refx))
    # # @test rel_err(dsol.q, refx) < 4E-6
    #
    # dint = Integrator(vdae, getTableauVSPARKGLRKLobIIIAC(2), Δt)
    # dsol = integrate(dint, nt)
    #
    # println(rel_err(dsol.q, refx))
    # # @test rel_err(dsol.q, refx) < 1E-11
    #
    # dint = Integrator(vdae, getTableauVSPARKGLRKLobIIIAC(3), Δt)
    # dsol = integrate(dint, nt)
    #
    # println(rel_err(dsol.q, refx))
    # # @test rel_err(dsol.q, refx) < 1E-11
    #
    # dint = Integrator(vdae, getTableauVSPARKGLRKLobIIIAD(1), Δt)
    # dsol = integrate(dint, nt)
    #
    # println(rel_err(dsol.q, refx))
    # # @test rel_err(dsol.q, refx) < 4E-6
    #
    # dint = Integrator(vdae, getTableauVSPARKGLRKLobIIIAD(2), Δt)
    # dsol = integrate(dint, nt)
    #
    # println(rel_err(dsol.q, refx))
    # # @test rel_err(dsol.q, refx) < 1E-11
    #
    # dint = Integrator(vdae, getTableauVSPARKGLRKLobIIIAD(3), Δt)
    # dsol = integrate(dint, nt)
    #
    # println(rel_err(dsol.q, refx))
    # # @test rel_err(dsol.q, refx) < 1E-11
    #
    # dint = Integrator(vdae, getTableauVSPARKGLRKLobIIIAE(1), Δt)
    # dsol = integrate(dint, nt)
    #
    # println(rel_err(dsol.q, refx))
    # # @test rel_err(dsol.q, refx) < 4E-6
    #
    # dint = Integrator(vdae, getTableauVSPARKGLRKLobIIIAE(2), Δt)
    # dsol = integrate(dint, nt)
    #
    # println(rel_err(dsol.q, refx))
    # # @test rel_err(dsol.q, refx) < 1E-11
    #
    # dint = Integrator(vdae, getTableauVSPARKGLRKLobIIIAE(3), Δt)
    # dsol = integrate(dint, nt)
    #
    # println(rel_err(dsol.q, refx))
    # # @test rel_err(dsol.q, refx) < 1E-11
end


# @testset "$(rpad("HPARK integrators",80))" begin
#
#     dint = Integrator(pdae, getTableauHPARK(:hpark_glrk1, getCoefficientsGLRK(1), getCoefficientsGLRK(1)), Δt)
#     dsol = integrate(dint, nt)
#
#     # println(rel_err(dsol.q, refx))
#     @test rel_err(dsol.q, refx) < 2E-6
#
#     dint = Integrator(pdae, getTableauHPARK(:hpark_glrk2, getCoefficientsGLRK(2), getCoefficientsGLRK(2)), Δt)
#     dsol = integrate(dint, nt)
#
#     # println(rel_err(dsol.q, refx))
#     @test rel_err(dsol.q, refx) < 8E-7
#
#     dint = Integrator(pdae, getTableauHPARKGLRK(1), Δt)
#     dsol = integrate(dint, nt)
#
#     # println(rel_err(dsol.q, refx))
#     @test rel_err(dsol.q, refx) < 2E-6
#
#     dint = Integrator(pdae, getTableauHPARKGLRK(2), Δt)
#     dsol = integrate(dint, nt)
#
#     # println(rel_err(dsol.q, refx))
#     @test rel_err(dsol.q, refx) < 8E-7
#
#     dint = Integrator(pdae, getTableauHPARKLobIIIAIIIB2(), Δt)
#     # dsol = integrate(dint, nt)
#     # TODO
#     # println(rel_err(dsol.q, refx))
#     # @test rel_err(dsol.q, refx) < 2E-2
#
#     dint = Integrator(pdae, getTableauHPARKLobIIIAIIIB3(), Δt)
#     # dsol = integrate(dint, nt)
#     # TODO
#     # println(rel_err(dsol.q, refx))
#     # @test rel_err(dsol.q, refx) < 8E-2
#
# end
#
#
# @testset "$(rpad("HSPARK integrators",80))" begin
#     ### HSPARK Integrators ###
#
#     dint = IntegratorHSPARK(pdae, getTableauSPARKGLRK(1), Δt)
#     dsol = integrate(dint, nt)
#
#     # println(rel_err(dsol.q, refx))
#     @test rel_err(dsol.q, refx) < 1E-6
#
#
#     ### HSPARKprimary Integrators ###
#
#     dint = Integrator(pdae, getTableauHSPARKGLRKpSymmetric(1), Δt)
#     dsol = integrate(dint, nt)
#
#     # println(rel_err(dsol.q, refx))
#     @test rel_err(dsol.q, refx) < 2E-1
#
#     dint = Integrator(pdae, getTableauHSPARKGLRKpSymmetric(2), Δt)
#     dsol = integrate(dint, nt)
#
#     # println(rel_err(dsol.q, refx))
#     @test rel_err(dsol.q, refx) < 4E-6
#
#     dint = Integrator(pdae, getTableauHSPARKLobIIIAIIIB2pSymmetric(), Δt)
#     dsol = integrate(dint, nt)
#
#     # println(rel_err(dsol.q, refx))
#     @test rel_err(dsol.q, refx) < 2E-1
#
#     dint = Integrator(pdae, getTableauHSPARKLobIIIAIIIB3pSymmetric(), Δt)
#     dsol = integrate(dint, nt)
#
#     # println(rel_err(dsol.q, refx))
#     @test rel_err(dsol.q, refx) < 4E-6
# end
