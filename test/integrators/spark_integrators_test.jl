module SPARKIntegratorsTest

    export test_spark_integrators

    using GeometricIntegrators
    using Test

    set_config(:nls_solver, NewtonSolver)
    set_config(:jacobian_autodiff, false)

    using ..OscillatorTest
    using ..OscillatorTest: Δt, nt, refx

    include("utils.jl")

    idae = oscillator_idae()
    pdae = oscillator_pdae()


    function test_spark_integrators()
        ### VPARK Integrator ###

        dint = Integrator(idae, getTableauSymplecticProjection(:pglrk1ps, getCoefficientsGLRK(1), getCoefficientsGLRK(1)), Δt)
        dsol = integrate(dint, nt)

        @test rel_err(dsol.q, refx) < 5E-4

        dint = Integrator(idae, getTableauSymplecticProjection(:pglrk2ps, getCoefficientsGLRK(2), getCoefficientsGLRK(2)), Δt)
        dsol = integrate(dint, nt)

        @test rel_err(dsol.q, refx) < 5E-8

        dint = Integrator(idae, getTableauGLRKpSymplectic(1), Δt)
        dsol = integrate(dint, nt)

        @test rel_err(dsol.q, refx) < 5E-4

        dint = Integrator(idae, getTableauGLRKpSymplectic(2), Δt)
        dsol = integrate(dint, nt)

        @test rel_err(dsol.q, refx) < 5E-8

        dint = Integrator(idae, getTableauLobIIIAIIIB2pSymplectic(), Δt)
        dsol = integrate(dint, nt)

        @test rel_err(dsol.q, refx) < 8E-4

        dint = Integrator(idae, getTableauLobIIIAIIIB3pSymplectic(), Δt)
        dsol = integrate(dint, nt)

        @test rel_err(dsol.q, refx) < 2E-4


        ### VSPARK Integrators ###

        dint = Integrator(idae, getTableauSymmetricProjection(:pglrk2ps, getCoefficientsGLRK(1), getCoefficientsGLRK(1)), Δt)
        dsol = integrate(dint, nt)

        @test rel_err(dsol.q, refx) < 5E-4

        dint = Integrator(idae, getTableauSymmetricProjection(:pglrk2ps, getCoefficientsGLRK(2), getCoefficientsGLRK(2)), Δt)
        dsol = integrate(dint, nt)

        @test rel_err(dsol.q, refx) < 5E-8

        dint = Integrator(idae, getTableauGLRKpSymmetric(1), Δt)
        dsol = integrate(dint, nt)

        @test rel_err(dsol.q, refx) < 5E-4

        dint = Integrator(idae, getTableauGLRKpSymmetric(2), Δt)
        dsol = integrate(dint, nt)

        @test rel_err(dsol.q, refx) < 5E-8

        dint = Integrator(idae, getTableauLobIIIAIIIB2pSymmetric(), Δt)
        dsol = integrate(dint, nt)

        @test rel_err(dsol.q, refx) < 8E-3

        # dint = Integrator(idae, getTableauLobIIIAIIIB3pSymmetric(), Δt)
        # dsol = integrate(dint, nt)
        #
        # println(rel_err(dsol.q, refx))
        # @test rel_err(dsol.q, refx) < 1E-2


        ### HPARK Integrator ###

        dint = Integrator(pdae, getTableauSymplecticProjection(:pglrk1ps, getCoefficientsGLRK(1), getCoefficientsGLRK(1)), Δt)
        dsol = integrate(dint, nt)

        @test rel_err(dsol.q, refx) < 5E-4

        dint = Integrator(pdae, getTableauSymplecticProjection(:pglrk2ps, getCoefficientsGLRK(2), getCoefficientsGLRK(2)), Δt)
        dsol = integrate(dint, nt)

        @test rel_err(dsol.q, refx) < 5E-8

        dint = Integrator(pdae, getTableauGLRKpSymplectic(1), Δt)
        dsol = integrate(dint, nt)

        @test rel_err(dsol.q, refx) < 5E-4

        dint = Integrator(pdae, getTableauGLRKpSymplectic(2), Δt)
        dsol = integrate(dint, nt)

        @test rel_err(dsol.q, refx) < 5E-8

        dint = Integrator(pdae, getTableauLobIIIAIIIB2pSymplectic(), Δt)
        dsol = integrate(dint, nt)

        @test rel_err(dsol.q, refx) < 2E-2

        dint = Integrator(pdae, getTableauLobIIIAIIIB3pSymplectic(), Δt)
        dsol = integrate(dint, nt)

        @test rel_err(dsol.q, refx) < 2E-4
    end
end

using .SPARKIntegratorsTest
test_spark_integrators()
