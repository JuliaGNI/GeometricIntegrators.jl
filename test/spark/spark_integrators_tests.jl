using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricProblems.LotkaVolterra2d
using RungeKutta
using Test

const t₀ = 0.0
const q₀ = [1.0, 1.0]
const params = (a₁=1.0, a₂=1.0, b₁=-1.0, b₂=-2.0)

const Δt = 0.01
const nt = 10
const tspan = (t₀, Δt * nt)

ode = odeproblem(q₀; tspan=tspan, tstep=Δt, parameters=params)
hdae = hdaeproblem(q₀; tspan=tspan, tstep=Δt, parameters=params)
idae = idaeproblem(q₀; tspan=tspan, tstep=Δt, parameters=params)
pdae = pdaeproblem(q₀; tspan=tspan, tstep=Δt, parameters=params)
ldae = ldaeproblem(q₀; tspan=tspan, tstep=Δt, parameters=params)
ldae_slrk = ldaeproblem_slrk(q₀; tspan=tspan, tstep=Δt, parameters=params)

ref = integrate(ode, Gauss(8))


@testset "$(rpad("SLRK integrators",80))" begin

    sol = integrate(ldae_slrk, SLRKLobattoIIIAB(2))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6

    sol = integrate(ldae_slrk, SLRKLobattoIIIAB(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    sol = integrate(ldae_slrk, SLRKLobattoIIIAB(4))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15


    sol = integrate(ldae_slrk, SLRKLobattoIIIBA(2))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6

    sol = integrate(ldae_slrk, SLRKLobattoIIIBA(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    sol = integrate(ldae_slrk, SLRKLobattoIIIBA(4))
    @test relative_maximum_error(sol.q, ref.q) < 1E-15


    sol = integrate(ldae_slrk, SLRKLobattoIIICC̄(2))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6

    sol = integrate(ldae_slrk, SLRKLobattoIIICC̄(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    sol = integrate(ldae_slrk, SLRKLobattoIIICC̄(4))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15


    sol = integrate(ldae_slrk, SLRKLobattoIIIC̄C(2))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6

    sol = integrate(ldae_slrk, SLRKLobattoIIIC̄C(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    sol = integrate(ldae_slrk, SLRKLobattoIIIC̄C(4))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15


    sol = integrate(ldae_slrk, SLRKLobattoIIID(2))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6

    sol = integrate(ldae_slrk, SLRKLobattoIIID(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    sol = integrate(ldae_slrk, SLRKLobattoIIID(4))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15


    sol = integrate(ldae_slrk, SLRKLobattoIIIE(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    sol = integrate(ldae_slrk, SLRKLobattoIIIE(3))
    @test relative_maximum_error(sol.q, ref.q) < 1E-11

    sol = integrate(ldae_slrk, SLRKLobattoIIIE(4))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15

end


@testset "$(rpad("SPARK integrators",80))" begin

    sol = integrate(idae, SPARKGLRK(1))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    sol = integrate(idae, SPARKGLRK(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-11


    # TODO: Check Errors !!!  GLVPRK should do much better ! (maybe problem with R∞?)

    sol = integrate(idae, SPARKGLVPRK(1))
    # println(relative_maximum_error(sol.q, ref.q))
    @test relative_maximum_error(sol.q, ref.q) < 2E-6

    sol = integrate(idae, SPARKGLVPRK(2))
    # println(relative_maximum_error(sol.q, ref.q))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6


    sol = integrate(idae, SPARKLobABC(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    sol = integrate(idae, SPARKLobABC(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    sol = integrate(idae, SPARKLobABC(4))
    @test relative_maximum_error(sol.q, ref.q) < 1E-15


    sol = integrate(idae, SPARKLobABD(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    sol = integrate(idae, SPARKLobABD(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    sol = integrate(idae, SPARKLobABD(4))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15


    # TODO: Check why these don't work !!!

    # sol = integrate(idae, SPARKLobattoIIIAIIIB(2))
    # println(relative_maximum_error(sol.q, ref.q))
    # @test relative_maximum_error(sol.q, ref.q) < 1E-6

    # sol = integrate(idae, SPARKLobattoIIIAIIIB(3))
    # println(relative_maximum_error(sol.q, ref.q))
    # @test relative_maximum_error(sol.q, ref.q) < 1E-6

    # sol = integrate(idae, SPARKLobattoIIIAIIIB(4))
    # println(relative_maximum_error(sol.q, ref.q))
    # @test relative_maximum_error(sol.q, ref.q) < 2E-10


    # TODO: Check why these don't work !!!

    # sol = integrate(idae, SPARKLobattoIIIBIIIA(2))
    # println(relative_maximum_error(sol.q, ref.q))
    # @test relative_maximum_error(sol.q, ref.q) < 1E-6

    # sol = integrate(idae, SPARKLobattoIIIBIIIA(3))
    # println(relative_maximum_error(sol.q, ref.q))
    # @test relative_maximum_error(sol.q, ref.q) < 1E-6

    # sol = integrate(idae, SPARKLobattoIIIBIIIA(4))
    # println(relative_maximum_error(sol.q, ref.q))
    # @test relative_maximum_error(sol.q, ref.q) < 2E-10


    # TODO: Check if the following integrators show the correct order of convergence !

    # sol = integrate(idae, SPARKGLRKLobattoIIIAIIIB(1))
    # # println(relative_maximum_error(sol.q, ref.q))
    # @test relative_maximum_error(sol.q, ref.q) < 5E-4

    # sol = integrate(idae, SPARKGLRKLobattoIIIAIIIB(2))
    # # println(relative_maximum_error(sol.q, ref.q))
    # @test relative_maximum_error(sol.q, ref.q) < 3E-4

    # sol = integrate(idae, SPARKGLRKLobattoIIIAIIIB(3))
    # # println(relative_maximum_error(sol.q, ref.q))
    # @test relative_maximum_error(sol.q, ref.q) < 2E-4


    # sol = integrate(idae, SPARKGLRKLobattoIIIBIIIA(1))
    # # println(relative_maximum_error(sol.q, ref.q))
    # @test relative_maximum_error(sol.q, ref.q) < 5E-4

    # sol = integrate(idae, SPARKGLRKLobattoIIIBIIIA(2))
    # # println(relative_maximum_error(sol.q, ref.q))
    # @test relative_maximum_error(sol.q, ref.q) < 3E-4

    # sol = integrate(idae, SPARKGLRKLobattoIIIBIIIA(3))
    # # println(relative_maximum_error(sol.q, ref.q))
    # @test relative_maximum_error(sol.q, ref.q) < 2E-4

end


@testset "$(rpad("VPARK integrators",80))" begin

    sol = integrate(idae, TableauSymplecticProjection(:glrk1ps, TableauGauss(1), TableauGauss(1)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    sol = integrate(idae, TableauSymplecticProjection(:glrk2ps, TableauGauss(2), TableauGauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-11

    sol = integrate(idae, TableauGausspSymplectic(1))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    sol = integrate(idae, TableauGausspSymplectic(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-11

    sol = integrate(idae, TableauLobattoIIIAIIIBpSymplectic(2))
    @test relative_maximum_error(sol.q, ref.q) < 2E-6

    sol = integrate(idae, TableauLobattoIIIAIIIBpSymplectic(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    sol = integrate(idae, TableauLobattoIIIAIIIBpSymplectic(4))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15

    sol = integrate(idae, TableauLobattoIIIBIIIApSymplectic(2))
    @test relative_maximum_error(sol.q, ref.q) < 2E-6

    sol = integrate(idae, TableauLobattoIIIBIIIApSymplectic(3))
    @test relative_maximum_error(sol.q, ref.q) < 4E-11

    sol = integrate(idae, TableauLobattoIIIBIIIApSymplectic(4))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15

end


# TODO: Replace idae with ldae !!!

@testset "$(rpad("VSPARK integrators",80))" begin

    # TODO: Fix or understand why these are not working !!!

    # sol = integrate(idae, VSPARK(SPARKGLRK(1)))
    # @test relative_maximum_error(sol.q, ref.q) < 1E-6

    # sol = integrate(idae, VSPARK(SPARKGLRK(2)))
    # @test relative_maximum_error(sol.q, ref.q) < 1E-11


    # TODO: Fix or understand why these are not working !!!

    # sol = integrate(idae, VSPARK(SPARKLobABC(2)))
    # @test relative_maximum_error(sol.q, ref.q) < 1E-6

    # sol = integrate(idae, VSPARK(SPARKLobABC(3)))
    # @test relative_maximum_error(sol.q, ref.q) < 2E-11

    # sol = integrate(idae, VSPARK(SPARKLobABC(4)))
    # @test relative_maximum_error(sol.q, ref.q) < 1E-15


    # sol = integrate(idae, VSPARK(SPARKLobABD(2)))
    # @test relative_maximum_error(sol.q, ref.q) < 1E-6

    # sol = integrate(idae, VSPARK(SPARKLobABD(3)))
    # @test relative_maximum_error(sol.q, ref.q) < 2E-11

    # sol = integrate(idae, VSPARK(SPARKLobABD(4)))
    # @test relative_maximum_error(sol.q, ref.q) < 2E-15


    # TODO: Fix or understand why these are not working !!!

    # sol = integrate(idae, VSPARK(SPARKLobattoIIIAIIIB(2)))
    # @test relative_maximum_error(sol.q, ref.q) < 1E-6

    # sol = integrate(idae, VSPARK(SPARKLobattoIIIAIIIB(3)))
    # @test relative_maximum_error(sol.q, ref.q) < 1E-6

    # sol = integrate(idae, VSPARK(SPARKLobattoIIIAIIIB(4)))
    # @test relative_maximum_error(sol.q, ref.q) < 2E-10


    # sol = integrate(idae, VSPARK(SPARKLobattoIIIBIIIA(2)))
    # @test relative_maximum_error(sol.q, ref.q) < 1E-6

    # sol = integrate(idae, VSPARK(SPARKLobattoIIIBIIIA(3)))
    # @test relative_maximum_error(sol.q, ref.q) < 1E-6

    # sol = integrate(idae, VSPARK(SPARKLobattoIIIBIIIA(4)))
    # @test relative_maximum_error(sol.q, ref.q) < 2E-10


    sol = integrate(idae, VSPARK(SPARKGLRKLobattoIIIAIIIB(1)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    sol = integrate(idae, VSPARK(SPARKGLRKLobattoIIIAIIIB(2)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-11

    sol = integrate(idae, VSPARK(SPARKGLRKLobattoIIIAIIIB(3)))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15


    sol = integrate(idae, VSPARK(SPARKGLRKLobattoIIIBIIIA(1)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    sol = integrate(idae, VSPARK(SPARKGLRKLobattoIIIBIIIA(2)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-11

    sol = integrate(idae, VSPARK(SPARKGLRKLobattoIIIBIIIA(3)))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15

end


@testset "$(rpad("VSPARK integrators with projection on primary constraint",80))" begin

    ## VSPARKprimary Integrators ###

    sol = integrate(idae, TableauVSPARKGLRKpMidpoint(1))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    sol = integrate(idae, TableauVSPARKGLRKpMidpoint(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-11

    sol = integrate(idae, TableauVSPARKGLRKpSymplectic(1))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    sol = integrate(idae, TableauVSPARKGLRKpSymplectic(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-11

    sol = integrate(idae, TableauVSPARKGLRKpSymmetric(1))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    sol = integrate(idae, TableauVSPARKGLRKpSymmetric(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-11

    sol = integrate(idae, TableauVSPARKLobattoIIIAIIIBpSymmetric(2))
    @test relative_maximum_error(sol.q, ref.q) < 2E-6

    sol = integrate(idae, TableauVSPARKLobattoIIIAIIIBpSymmetric(3))
    @test relative_maximum_error(sol.q, ref.q) < 5E-11

    sol = integrate(idae, TableauVSPARKLobattoIIIAIIIBpSymmetric(4))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15

    sol = integrate(idae, TableauVSPARKLobattoIIIBIIIApSymmetric(2))
    @test relative_maximum_error(sol.q, ref.q) < 2E-6

    sol = integrate(idae, TableauVSPARKLobattoIIIBIIIApSymmetric(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-5
    # @test relative_maximum_error(sol.q, ref.q) < 5E-11
    # TODO: Check Errors !!!

    sol = integrate(idae, TableauVSPARKLobattoIIIBIIIApSymmetric(4))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15

end


@testset "$(rpad("VSPARK integrators with projection on secondary constraint",80))" begin

    ## VSPARKsecondary Integrators ###

    sol = integrate(ldae, TableauVSPARKLobattoIIIAB(2))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6

    sol = integrate(ldae, TableauVSPARKLobattoIIIAB(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    sol = integrate(ldae, TableauVSPARKLobattoIIIAB(4))
    @test relative_maximum_error(sol.q, ref.q) < 4E-15


    sol = integrate(ldae, TableauVSPARKLobattoIIIBA(2))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6

    sol = integrate(ldae, TableauVSPARKLobattoIIIBA(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    sol = integrate(ldae, TableauVSPARKLobattoIIIBA(4))
    @test relative_maximum_error(sol.q, ref.q) < 1E-15


    sol = integrate(ldae, TableauVSPARKLobattoIIICC̄(2))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6

    sol = integrate(ldae, TableauVSPARKLobattoIIICC̄(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    sol = integrate(ldae, TableauVSPARKLobattoIIICC̄(4))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15


    sol = integrate(ldae, TableauVSPARKLobattoIIIC̄C(2))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6

    sol = integrate(ldae, TableauVSPARKLobattoIIIC̄C(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    sol = integrate(ldae, TableauVSPARKLobattoIIIC̄C(4))
    @test relative_maximum_error(sol.q, ref.q) < 1E-15


    sol = integrate(ldae, TableauVSPARKLobattoIIID(2))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6

    sol = integrate(ldae, TableauVSPARKLobattoIIID(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    sol = integrate(ldae, TableauVSPARKLobattoIIID(4))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15


    sol = integrate(ldae, TableauVSPARKLobattoIIIE(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    sol = integrate(ldae, TableauVSPARKLobattoIIIE(3))
    @test relative_maximum_error(sol.q, ref.q) < 1E-11

    sol = integrate(ldae, TableauVSPARKLobattoIIIE(4))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15


    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIIAB(1))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6

    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIIAB(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-11

    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIIAB(3))
    @test relative_maximum_error(sol.q, ref.q) < 1E-15


    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIIBA(1))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6

    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIIBA(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-11

    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIIBA(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15


    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIICC̄(1))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6

    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIICC̄(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-11

    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIICC̄(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15


    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIIC̄C(1))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6

    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIIC̄C(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-11

    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIIC̄C(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15


    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIID(1))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6

    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIID(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-11

    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIID(3))
    @test relative_maximum_error(sol.q, ref.q) < 4E-15
    # @test relative_maximum_error(sol.q, ref.q) < 2E-15
    # TODO: Check errors!


    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIIE(1))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIIE(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-11

    sol = integrate(ldae, TableauVSPARKGLRKLobattoIIIE(3))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15

end


@testset "$(rpad("HPARK integrators",80))" begin

    sol = integrate(pdae, TableauHPARKGLRK(1))
    # println(relative_maximum_error(sol.q, ref.q))
    @test relative_maximum_error(sol.q, ref.q) < 2E-6

    sol = integrate(pdae, TableauHPARKGLRK(2))
    # println(relative_maximum_error(sol.q, ref.q))
    @test relative_maximum_error(sol.q, ref.q) < 8E-7

    # sol = integrate(pdae, TableauHPARKLobattoIIIAIIIB(2))
    # println(relative_maximum_error(sol.q, ref.q))
    # @test relative_maximum_error(sol.q, ref.q) < 2E-2
    # TODO: Check errors and large number of solver iterations !!!

    # sol = integrate(pdae, TableauHPARKLobattoIIIAIIIB(3))
    # println(relative_maximum_error(sol.q, ref.q))
    # @test relative_maximum_error(sol.q, ref.q) < 8E-2
    # TODO: Check errors and large number of solver iterations !!!

    sol = integrate(pdae, TableauHPARKLobattoIIIAIIIB(4))
    # println(relative_maximum_error(sol.q, ref.q))
    @test relative_maximum_error(sol.q, ref.q) < 2E-3
    # TODO: Check errors and large number of solver iterations !!!

    # sol = integrate(pdae, TableauHPARKLobattoIIIBIIIA(2))
    # println(relative_maximum_error(sol.q, ref.q))
    # @test relative_maximum_error(sol.q, ref.q) < 2E-2
    # TODO: Check errors and large number of solver iterations !!!

    # sol = integrate(pdae, TableauHPARKLobattoIIIBIIIA(3))
    # println(relative_maximum_error(sol.q, ref.q))
    # @test relative_maximum_error(sol.q, ref.q) < 4E-3
    # TODO: Check errors and large number of solver iterations !!!

    sol = integrate(pdae, TableauHPARKLobattoIIIBIIIA(4))
    # println(relative_maximum_error(sol.q, ref.q))
    @test relative_maximum_error(sol.q, ref.q) < 4E-2
    # TODO: Check errors and large number of solver iterations !!!

end


# TODO: Replace pdae with hdae !!!

@testset "$(rpad("HSPARK integrators",80))" begin

    # println("HSPARK(SPARKGLRK(1))")
    sol = integrate(pdae, HSPARK(SPARKGLRK(1)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    # println("HSPARK(SPARKGLRK(2))")
    sol = integrate(pdae, HSPARK(SPARKGLRK(2)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-11


    # println("HSPARK(SPARKLobABC(2))")
    sol = integrate(pdae, HSPARK(SPARKLobABC(2)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    # println("HSPARK(SPARKLobABC(3))")
    sol = integrate(pdae, HSPARK(SPARKLobABC(3)))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    # println("HSPARK(SPARKLobABC(4))")
    sol = integrate(pdae, HSPARK(SPARKLobABC(4)))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15


    # println("HSPARK(SPARKLobABD(2))")
    sol = integrate(pdae, HSPARK(SPARKLobABD(2)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-6

    # println("HSPARK(SPARKLobABD(3))")
    sol = integrate(pdae, HSPARK(SPARKLobABD(3)))
    @test relative_maximum_error(sol.q, ref.q) < 2E-11

    # println("HSPARK(SPARKLobABD(4))")
    sol = integrate(pdae, HSPARK(SPARKLobABD(4)))
    @test relative_maximum_error(sol.q, ref.q) < 2E-15


    # TODO: Check why these don't work properly !!!

    # sol = integrate(pdae, HSPARK(SPARKLobattoIIIAIIIB(2)))
    # println(relative_maximum_error(sol.q, ref.q))
    # @test relative_maximum_error(sol.q, ref.q) < 1E-6

    # sol = integrate(pdae, HSPARK(SPARKLobattoIIIAIIIB(3)))
    # println(relative_maximum_error(sol.q, ref.q))
    # @test relative_maximum_error(sol.q, ref.q) < 1E-6

    # sol = integrate(pdae, HSPARK(SPARKLobattoIIIAIIIB(4)))
    # println(relative_maximum_error(sol.q, ref.q))
    # @test relative_maximum_error(sol.q, ref.q) < 2E-10


    # TODO: Check why these don't work properly !!!

    # sol = integrate(pdae, HSPARK(SPARKLobattoIIIBIIIA(2)))
    # println(relative_maximum_error(sol.q, ref.q))
    # @test relative_maximum_error(sol.q, ref.q) < 1E-6

    # sol = integrate(pdae, HSPARK(SPARKLobattoIIIBIIIA(3)))
    # println(relative_maximum_error(sol.q, ref.q))
    # @test relative_maximum_error(sol.q, ref.q) < 1E-6

    # sol = integrate(pdae, HSPARK(SPARKLobattoIIIBIIIA(4)))
    # println(relative_maximum_error(sol.q, ref.q))
    # @test relative_maximum_error(sol.q, ref.q) < 2E-10


    # TODO: Check if the following integrators show the correct order of convergence !

    # sol = integrate(pdae, HSPARK(SPARKGLRKLobattoIIIAIIIB(1)))
    # println(relative_maximum_error(sol.q, ref.q))
    # @test relative_maximum_error(sol.q, ref.q) < 5E-4

    # sol = integrate(pdae, HSPARK(SPARKGLRKLobattoIIIAIIIB(2)))
    # println(relative_maximum_error(sol.q, ref.q))
    # @test relative_maximum_error(sol.q, ref.q) < 3E-4

    # sol = integrate(pdae, HSPARK(SPARKGLRKLobattoIIIAIIIB(3)))
    # println(relative_maximum_error(sol.q, ref.q))
    # @test relative_maximum_error(sol.q, ref.q) < 2E-4


    # sol = integrate(pdae, HSPARK(SPARKGLRKLobattoIIIBIIIA(1)))
    # println(relative_maximum_error(sol.q, ref.q))
    # @test relative_maximum_error(sol.q, ref.q) < 5E-4

    # sol = integrate(pdae, HSPARK(SPARKGLRKLobattoIIIBIIIA(2)))
    # println(relative_maximum_error(sol.q, ref.q))
    # @test relative_maximum_error(sol.q, ref.q) < 3E-4

    # sol = integrate(pdae, HSPARK(SPARKGLRKLobattoIIIBIIIA(3)))
    # println(relative_maximum_error(sol.q, ref.q))
    # @test relative_maximum_error(sol.q, ref.q) < 2E-4

end


@testset "$(rpad("HSPARK integrators with projection on primary constraint",80))" begin

    ### HSPARKprimary Integrators ###

    sol = integrate(pdae, TableauHSPARKGLRKpSymmetric(1))
    # println(relative_maximum_error(sol.q, ref.q))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6
    # TODO: Check Errors !!!

    sol = integrate(pdae, TableauHSPARKGLRKpSymmetric(2))
    # println(relative_maximum_error(sol.q, ref.q))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6
    # TODO: Check Errors !!!

    sol = integrate(pdae, TableauHSPARKLobattoIIIAIIIBpSymmetric(2))
    # println(relative_maximum_error(sol.q, ref.q))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6
    # TODO: Check Errors !!!

    sol = integrate(pdae, TableauHSPARKLobattoIIIAIIIBpSymmetric(3))
    # println(relative_maximum_error(sol.q, ref.q))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6
    # TODO: Check Errors !!!

    sol = integrate(pdae, TableauHSPARKLobattoIIIAIIIBpSymmetric(4))
    # println(relative_maximum_error(sol.q, ref.q))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6
    # TODO: Check Errors !!!

    sol = integrate(pdae, TableauHSPARKLobattoIIIBIIIApSymmetric(2))
    # println(relative_maximum_error(sol.q, ref.q))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6
    # TODO: Check Errors !!!

    sol = integrate(pdae, TableauHSPARKLobattoIIIBIIIApSymmetric(3))
    # println(relative_maximum_error(sol.q, ref.q))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6
    # TODO: Check Errors !!!

    sol = integrate(pdae, TableauHSPARKLobattoIIIBIIIApSymmetric(4))
    # println(relative_maximum_error(sol.q, ref.q))
    @test relative_maximum_error(sol.q, ref.q) < 4E-6
    # TODO: Check Errors !!!

end


# @testset "$(rpad("HSPARK integrators with projection on secondary constraint",80))" begin

#     ### HSPARKsecondary Integrators ###

#     sol = integrate(hdae, TableauHSPARKLobattoIIIAB(2))
#     # println(relative_maximum_error(sol.q, ref.q))
#     @test relative_maximum_error(sol.q, ref.q) < 1E-6

#     sol = integrate(hdae, TableauHSPARKLobattoIIIAB(3))
#     # println(relative_maximum_error(sol.q, ref.q))
#     # @test relative_maximum_error(sol.q, ref.q) < 2E-11
#     # TODO: Check Errors !!!

#     sol = integrate(hdae, TableauHSPARKLobattoIIIAB(4))
#     # println(relative_maximum_error(sol.q, ref.q))
#     # @test relative_maximum_error(sol.q, ref.q) < 1E-15
#     # TODO: Check Errors !!!


#     sol = integrate(hdae, TableauHSPARKLobattoIIIBA(2))
#     # println(relative_maximum_error(sol.q, ref.q))
#     @test relative_maximum_error(sol.q, ref.q) < 1E-6

#     sol = integrate(hdae, TableauHSPARKLobattoIIIBA(3))
#     # println(relative_maximum_error(sol.q, ref.q))
#     # @test relative_maximum_error(sol.q, ref.q) < 2E-11
#     # TODO: Check Errors !!!

#     sol = integrate(hdae, TableauHSPARKLobattoIIIBA(4))
#     # println(relative_maximum_error(sol.q, ref.q))
#     # @test relative_maximum_error(sol.q, ref.q) < 1E-15
#     # TODO: Check Errors !!!


#     sol = integrate(hdae, TableauHSPARKLobattoIIID(2))
#     # println(relative_maximum_error(sol.q, ref.q))
#     @test relative_maximum_error(sol.q, ref.q) < 1E-6

#     sol = integrate(hdae, TableauHSPARKLobattoIIID(3))
#     # println(relative_maximum_error(sol.q, ref.q))
#     # @test relative_maximum_error(sol.q, ref.q) < 2E-11
#     # TODO: Check Errors !!!

#     sol = integrate(hdae, TableauHSPARKLobattoIIID(4))
#     # println(relative_maximum_error(sol.q, ref.q))
#     # @test relative_maximum_error(sol.q, ref.q) < 1E-15
#     # TODO: Check Errors !!!


#     sol = integrate(hdae, TableauHSPARKLobattoIIIE(2))
#     # println(relative_maximum_error(sol.q, ref.q))
#     @test relative_maximum_error(sol.q, ref.q) < 1E-6

#     sol = integrate(hdae, TableauHSPARKLobattoIIIE(3))
#     # println(relative_maximum_error(sol.q, ref.q))
#     # @test relative_maximum_error(sol.q, ref.q) < 2E-11
#     # TODO: Check Errors !!!

#     sol = integrate(hdae, TableauHSPARKLobattoIIIE(4))
#     # println(relative_maximum_error(sol.q, ref.q))
#     # @test relative_maximum_error(sol.q, ref.q) < 1E-15
#     # TODO: Check Errors !!!


#     # TODO: Fix the following !!!

#     # sol = integrate(hdae, TableauHSPARKGLRKLobattoIIIAB(2))
#     # println(relative_maximum_error(sol.q, ref.q))
#     # @test relative_maximum_error(sol.q, ref.q) < 4E-6

#     # sol = integrate(hdae, TableauHSPARKGLRKLobattoIIIAB(3))
#     # println(relative_maximum_error(sol.q, ref.q))
#     # @test relative_maximum_error(sol.q, ref.q) < 2E-11

#     # sol = integrate(hdae, TableauHSPARKGLRKLobattoIIIAB(4))
#     # println(relative_maximum_error(sol.q, ref.q))
#     # @test relative_maximum_error(sol.q, ref.q) < 1E-15


#     # sol = integrate(hdae, TableauHSPARKGLRKLobattoIIIBA(2))
#     # println(relative_maximum_error(sol.q, ref.q))
#     # @test relative_maximum_error(sol.q, ref.q) < 4E-6

#     # sol = integrate(hdae, TableauHSPARKGLRKLobattoIIIBA(3))
#     # println(relative_maximum_error(sol.q, ref.q))
#     # @test relative_maximum_error(sol.q, ref.q) < 2E-11

#     # sol = integrate(hdae, TableauHSPARKGLRKLobattoIIIBA(4))
#     # println(relative_maximum_error(sol.q, ref.q))
#     # @test relative_maximum_error(sol.q, ref.q) < 1E-15


#     # sol = integrate(hdae, TableauHSPARKGLRKLobattoIIID(2))
#     # println(relative_maximum_error(sol.q, ref.q))
#     # @test relative_maximum_error(sol.q, ref.q) < 4E-6

#     # sol = integrate(hdae, TableauHSPARKGLRKLobattoIIID(3))
#     # println(relative_maximum_error(sol.q, ref.q))
#     # @test relative_maximum_error(sol.q, ref.q) < 2E-11

#     # sol = integrate(hdae, TableauHSPARKGLRKLobattoIIID(4))
#     # println(relative_maximum_error(sol.q, ref.q))
#     # @test relative_maximum_error(sol.q, ref.q) < 1E-15


#     # sol = integrate(hdae, TableauHSPARKGLRKLobattoIIIE(2))
#     # println(relative_maximum_error(sol.q, ref.q))
#     # @test relative_maximum_error(sol.q, ref.q) < 4E-6

#     # sol = integrate(hdae, TableauHSPARKGLRKLobattoIIIE(3))
#     # println(relative_maximum_error(sol.q, ref.q))
#     # @test relative_maximum_error(sol.q, ref.q) < 2E-11

#     # sol = integrate(hdae, TableauHSPARKGLRKLobattoIIIE(4))
#     # println(relative_maximum_error(sol.q, ref.q))
#     # @test relative_maximum_error(sol.q, ref.q) < 1E-15

# end
