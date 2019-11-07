module SPARKIntegratorsTest

    export test_spark_integrators

    using GeometricIntegrators
    using GeometricIntegrators.Utils
    using Test

    set_config(:nls_atol, 8eps())
    set_config(:nls_rtol, 2eps())
    # set_config(:nls_solver, NewtonSolver)
    # set_config(:jacobian_autodiff, false)

    using ..LotkaVolterraTest
    using ..LotkaVolterraTest: Δt, nt

    idae = lotka_volterra_2d_idae()
    pdae = lotka_volterra_2d_pdae()

    int = IntegratorFIRK(lotka_volterra_2d_ode(), getTableauGLRK(8), Δt)
    sol = integrate(int, nt)
    refx = sol.q[:,end]


    function test_spark_integrators()

        ### VPARK Integrator ###

        dint = Integrator(idae, getTableauSymplecticProjection(:pglrk1ps, getCoefficientsGLRK(1), getCoefficientsGLRK(1)), Δt)
        dsol = integrate(dint, nt)

        # println(rel_err(dsol.q, refx))
        @test rel_err(dsol.q, refx) < 1E-6

        dint = Integrator(idae, getTableauSymplecticProjection(:pglrk2ps, getCoefficientsGLRK(2), getCoefficientsGLRK(2)), Δt)
        dsol = integrate(dint, nt)

        # println(rel_err(dsol.q, refx))
        @test rel_err(dsol.q, refx) < 1E-11

        dint = Integrator(idae, getTableauGLRKpSymplectic(1), Δt)
        dsol = integrate(dint, nt)

        # println(rel_err(dsol.q, refx))
        @test rel_err(dsol.q, refx) < 1E-6

        dint = Integrator(idae, getTableauGLRKpSymplectic(2), Δt)
        dsol = integrate(dint, nt)

        # println(rel_err(dsol.q, refx))
        @test rel_err(dsol.q, refx) < 1E-11

        dint = Integrator(idae, getTableauLobIIIAIIIB2pSymplectic(), Δt)
        dsol = integrate(dint, nt)

        # println(rel_err(dsol.q, refx))
        @test rel_err(dsol.q, refx) < 2E-6

        dint = Integrator(idae, getTableauLobIIIAIIIB3pSymplectic(), Δt)
        dsol = integrate(dint, nt)

        # println(rel_err(dsol.q, refx))
        @test rel_err(dsol.q, refx) < 8E-5


        ### VSPARK Integrators ###

        dint = Integrator(idae, getTableauGLRKpMidpoint(1), Δt)
        dsol = integrate(dint, nt)

        # println(rel_err(dsol.q, refx))
        @test rel_err(dsol.q, refx) < 1E-6

        dint = Integrator(idae, getTableauGLRKpMidpoint(2), Δt)
        dsol = integrate(dint, nt)

        # println(rel_err(dsol.q, refx))
        @test rel_err(dsol.q, refx) < 1E-11

        # dint = Integrator(idae, getTableauVSPARKGLRK(1), Δt)
        # dsol = integrate(dint, nt)
        #
        # println(rel_err(dsol.q, refx))
        # # @test rel_err(dsol.q, refx) < 1E-6
        #
        # dint = Integrator(idae, getTableauVSPARKGLRK(2), Δt)
        # dsol = integrate(dint, nt)
        #
        # println(rel_err(dsol.q, refx))
        # # @test rel_err(dsol.q, refx) < 1E-11


        ### VSPARKprimary Integrators ###

        dint = Integrator(idae, getTableauSymmetricProjection(:pglrk1dps, getCoefficientsGLRK(1), getCoefficientsGLRK(1); R∞=-1), Δt)
        dsol = integrate(dint, nt)

        # println(rel_err(dsol.q, refx))
        @test rel_err(dsol.q, refx) < 1E-6

        dint = Integrator(idae, getTableauSymmetricProjection(:pglrk2dps, getCoefficientsGLRK(2), getCoefficientsGLRK(2); R∞=+1), Δt)
        dsol = integrate(dint, nt)

        # println(rel_err(dsol.q, refx))
        @test rel_err(dsol.q, refx) < 1E-11

        dint = Integrator(idae, getTableauGLRKpSymmetric(1), Δt)
        dsol = integrate(dint, nt)

        # println(rel_err(dsol.q, refx))
        @test rel_err(dsol.q, refx) < 1E-6

        dint = Integrator(idae, getTableauGLRKpSymmetric(2), Δt)
        dsol = integrate(dint, nt)

        # println(rel_err(dsol.q, refx))
        @test rel_err(dsol.q, refx) < 1E-11

        dint = Integrator(idae, getTableauVSPARKGLRKpSymplectic(1), Δt)
        dsol = integrate(dint, nt)

        # println(rel_err(dsol.q, refx))
        @test rel_err(dsol.q, refx) < 1E-6

        dint = Integrator(idae, getTableauVSPARKGLRKpSymplectic(2), Δt)
        dsol = integrate(dint, nt)

        # println(rel_err(dsol.q, refx))
        @test rel_err(dsol.q, refx) < 1E-11

        dint = IntegratorVSPARKprimary(idae, getTableauLobIIIAIIIB2pSymmetric(), Δt)
        dsol = integrate(dint, nt)

        # println(rel_err(dsol.q, refx))
        @test rel_err(dsol.q, refx) < 2E-6

        dint = IntegratorVSPARKprimary(idae, getTableauLobIIIAIIIB3pSymmetric(), Δt)
        dsol = integrate(dint, nt)

        # println(rel_err(dsol.q, refx))
        @test rel_err(dsol.q, refx) < 8E-5


        ### HPARK Integrator ###

        dint = Integrator(pdae, getTableauHPARK(:hpark_glrk1, getCoefficientsGLRK(1), getCoefficientsGLRK(1)), Δt)
        dsol = integrate(dint, nt)

        # println(rel_err(dsol.q, refx))
        @test rel_err(dsol.q, refx) < 2E-6

        dint = Integrator(pdae, getTableauHPARK(:hpark_glrk2, getCoefficientsGLRK(2), getCoefficientsGLRK(2)), Δt)
        dsol = integrate(dint, nt)

        # println(rel_err(dsol.q, refx))
        @test rel_err(dsol.q, refx) < 8E-7

        dint = Integrator(pdae, getTableauHPARKGLRK(1), Δt)
        dsol = integrate(dint, nt)

        # println(rel_err(dsol.q, refx))
        @test rel_err(dsol.q, refx) < 2E-6

        dint = Integrator(pdae, getTableauHPARKGLRK(2), Δt)
        dsol = integrate(dint, nt)

        # println(rel_err(dsol.q, refx))
        @test rel_err(dsol.q, refx) < 8E-7

        # dint = Integrator(pdae, getTableauHPARKLobIIIAIIIB2(), Δt)
        # dsol = integrate(dint, nt)
        # TODO
        # println(rel_err(dsol.q, refx))
        # @test rel_err(dsol.q, refx) < 2E-2

        # dint = Integrator(pdae, getTableauHPARKLobIIIAIIIB3(), Δt)
        # dsol = integrate(dint, nt)
        # TODO
        # println(rel_err(dsol.q, refx))
        # @test rel_err(dsol.q, refx) < 2E-4


        ### HSPARK Integrators ###

        # set_config(:nls_nmax, 1)

        # dint = Integrator(pdae, getTableauSymmetricProjection(:pglrk2ps, getCoefficientsGLRK(1), getCoefficientsGLRK(1)), Δt)
        # dsol = integrate(dint, nt)
        #
        # println(rel_err(dsol.q, refx))
        # # @test rel_err(dsol.q, refx) < 5E-4
        #
        # dint = Integrator(pdae, getTableauSymmetricProjection(:pglrk2ps, getCoefficientsGLRK(2), getCoefficientsGLRK(2)), Δt)
        # dsol = integrate(dint, nt)
        #
        # println(rel_err(dsol.q, refx))
        # # @test rel_err(dsol.q, refx) < 5E-8
        #
        # dint = Integrator(pdae, getTableauGLRKpSymmetric(1), Δt)
        # dsol = integrate(dint, nt)
        #
        # println(rel_err(dsol.q, refx))
        # # @test rel_err(dsol.q, refx) < 5E-4
        #
        # dint = Integrator(pdae, getTableauGLRKpSymmetric(2), Δt)
        # dsol = integrate(dint, nt)
        #
        # println(rel_err(dsol.q, refx))
        # # @test rel_err(dsol.q, refx) < 5E-8

        # dint = Integrator(pdae, getTableauLobIIIAIIIB2pSymmetric(), Δt)
        # dsol = integrate(dint, nt)
        #
        # println(rel_err(dsol.q, refx))
        # @test rel_err(dsol.q, refx) < 8E-3
        #
        # dint = Integrator(pdae, getTableauLobIIIAIIIB3pSymmetric(), Δt)
        # dsol = integrate(dint, nt)
        #
        # println(rel_err(dsol.q, refx))
        # @test rel_err(dsol.q, refx) < 1E-2
    end
end

using .SPARKIntegratorsTest
test_spark_integrators()
