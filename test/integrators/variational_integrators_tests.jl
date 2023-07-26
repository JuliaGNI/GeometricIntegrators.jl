using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
using Test


iode = iodeproblem()
pref = exact_solution(podeproblem())


@testset "$(rpad("Vartiational Partitioned Runge-Kutta integrators",80))" begin

    sol = integrate(iode, VPRKGauss(1))
    # println(relative_maximum_error(sol.q, pref.q))
    # println(relative_maximum_error(sol.p, pref.p))
    # println()
    @test relative_maximum_error(sol.q, pref.q) < 4E-4
    # @test relative_maximum_error(sol.p, pref.p) < 4E-4

    sol = integrate(iode, VPRKGauss(2))
    # println(relative_maximum_error(sol.q, pref.q))
    # println(relative_maximum_error(sol.p, pref.p))
    # println()
    @test relative_maximum_error(sol.q, pref.q) < 4E-8
    # @test relative_maximum_error(sol.p, pref.p) < 4E-8

    sol = integrate(iode, VPRKGauss(3))
    # println(relative_maximum_error(sol.q, pref.q))
    # println(relative_maximum_error(sol.p, pref.p))
    # println()
    @test relative_maximum_error(sol.q, pref.q) < 8E-13
    # @test relative_maximum_error(sol.p, pref.p) < 8E-13

    sol = integrate(iode, VPRKLobattoIIIAIIIĀ(2))
    # println(relative_maximum_error(sol.q, pref.q))
    # println(relative_maximum_error(sol.p, pref.p))
    # println()
    @test relative_maximum_error(sol.q, pref.q) < 2E-4
    # @test relative_maximum_error(sol.p, pref.p) < 2E-4

    sol = integrate(iode, VPRKLobattoIIIAIIIĀ(3))
    # println(relative_maximum_error(sol.q, pref.q))
    # println(relative_maximum_error(sol.p, pref.p))
    # println()
    @test relative_maximum_error(sol.q, pref.q) < 8E-9
    # @test relative_maximum_error(sol.p, pref.p) < 8E-9

    sol = integrate(iode, VPRKLobattoIIIAIIIĀ(4))
    # println(relative_maximum_error(sol.q, pref.q))
    # println(relative_maximum_error(sol.p, pref.p))
    # println()
    @test relative_maximum_error(sol.q, pref.q) < 2E-13
    # @test relative_maximum_error(sol.p, pref.p) < 2E-13

    sol = integrate(iode, VPRKLobattoIIIBIIIB̄(2))
    # println(relative_maximum_error(sol.q, pref.q))
    # println(relative_maximum_error(sol.p, pref.p))
    # println()
    @test relative_maximum_error(sol.q, pref.q) < 4E-4
    # @test relative_maximum_error(sol.p, pref.p) < 4E-4

    sol = integrate(iode, VPRKLobattoIIIBIIIB̄(3))
    # println(relative_maximum_error(sol.q, pref.q))
    # println(relative_maximum_error(sol.p, pref.p))
    # println()
    @test relative_maximum_error(sol.q, pref.q) < 4E-8
    # @test relative_maximum_error(sol.p, pref.p) < 4E-8

    sol = integrate(iode, VPRKLobattoIIIBIIIB̄(4))
    # println(relative_maximum_error(sol.q, pref.q))
    # println(relative_maximum_error(sol.p, pref.p))
    # println()
    @test relative_maximum_error(sol.q, pref.q) < 8E-13
    # @test relative_maximum_error(sol.p, pref.p) < 8E-13

end
