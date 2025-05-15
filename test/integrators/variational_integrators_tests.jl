using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
using Test

using GeometricEquations: DELEProblem


lode = lodeproblem()
pref = exact_solution(podeproblem())


@testset "$(rpad("Discrete Euler-Lagrange integrators",80))" begin

    sol = integrate(deleproblem_midpoint(), DiscreteEulerLagrange())
    # println(relative_maximum_error(sol.q, pref.q))
    @test relative_maximum_error(sol.q, pref.q) < 4E-4

    dele = DELEProblem(lode, Midpoint())
    dsol = integrate(dele, DiscreteEulerLagrange())
    @test relative_maximum_error(dsol.q, pref.q) < 4E-4
    @test relative_maximum_error(dsol.q, sol.q) < 1E-14

    sol = integrate(deleproblem_trapezoidal(), DiscreteEulerLagrange())
    # println(relative_maximum_error(sol.q, pref.q))
    @test relative_maximum_error(sol.q, pref.q) < 4E-4

    dele = DELEProblem(lode, Trapezoidal())
    dsol = integrate(dele, DiscreteEulerLagrange())
    @test relative_maximum_error(dsol.q, pref.q) < 4E-4
    @test relative_maximum_error(dsol.q, sol.q) < 1E-14
end


@testset "$(rpad("Vartiational integrators",80))" begin

    sol = integrate(lode, PMVImidpoint())
    @test relative_maximum_error(sol.q, pref.q) < 4E-4

    sol = integrate(lode, PMVItrapezoidal())
    @test relative_maximum_error(sol.q, pref.q) < 4E-4

end


@testset "$(rpad("Vartiational Partitioned Runge-Kutta integrators",80))" begin

    sol = integrate(lode, VPRKGauss(1))
    # println(relative_maximum_error(sol.q, pref.q))
    # println(relative_maximum_error(sol.p, pref.p))
    # println()
    @test relative_maximum_error(sol.q, pref.q) < 4E-4
    # @test relative_maximum_error(sol.p, pref.p) < 4E-4

    ref = integrate(lode, PMVImidpoint())
    @test relative_maximum_error(sol.q, ref.q) < 8*eps()

    sol = integrate(lode, VPRKGauss(2))
    # println(relative_maximum_error(sol.q, pref.q))
    # println(relative_maximum_error(sol.p, pref.p))
    # println()
    @test relative_maximum_error(sol.q, pref.q) < 4E-8
    # @test relative_maximum_error(sol.p, pref.p) < 4E-8

    sol = integrate(lode, VPRKGauss(3))
    # println(relative_maximum_error(sol.q, pref.q))
    # println(relative_maximum_error(sol.p, pref.p))
    # println()
    @test relative_maximum_error(sol.q, pref.q) < 8E-13
    # @test relative_maximum_error(sol.p, pref.p) < 8E-13

    sol = integrate(lode, VPRKGauss(4))
    # println(relative_maximum_error(sol.q, pref.q))
    # println(relative_maximum_error(sol.p, pref.p))
    # println()
    @test relative_maximum_error(sol.q, pref.q) < 4E-16
    # @test relative_maximum_error(sol.q, pref.q) < 2E-16 # TODO: Reactivate!
    # @test relative_maximum_error(sol.p, pref.p) < 2E-16


    sol = integrate(lode, VPRKLobattoIIIAIIIĀ(2))
    # println(relative_maximum_error(sol.q, pref.q))
    # println(relative_maximum_error(sol.p, pref.p))
    # println()
    @test relative_maximum_error(sol.q, pref.q) < 2E-4
    # @test relative_maximum_error(sol.p, pref.p) < 2E-4

    ref = integrate(lode, PMVItrapezoidal())
    @test relative_maximum_error(sol.q, ref.q) < 8*eps()

    sol = integrate(lode, VPRKLobattoIIIAIIIĀ(3))
    # println(relative_maximum_error(sol.q, pref.q))
    # println(relative_maximum_error(sol.p, pref.p))
    # println()
    @test relative_maximum_error(sol.q, pref.q) < 8E-9
    # @test relative_maximum_error(sol.p, pref.p) < 8E-9

    sol = integrate(lode, VPRKLobattoIIIAIIIĀ(4))
    # println(relative_maximum_error(sol.q, pref.q))
    # println(relative_maximum_error(sol.p, pref.p))
    # println()
    @test relative_maximum_error(sol.q, pref.q) < 2E-13
    # @test relative_maximum_error(sol.p, pref.p) < 2E-13


    sol = integrate(lode, VPRKLobattoIIIBIIIB̄(2))
    # println(relative_maximum_error(sol.q, pref.q))
    # println(relative_maximum_error(sol.p, pref.p))
    # println()
    @test relative_maximum_error(sol.q, pref.q) < 4E-4
    # @test relative_maximum_error(sol.p, pref.p) < 4E-4

    sol = integrate(lode, VPRKLobattoIIIBIIIB̄(3))
    # println(relative_maximum_error(sol.q, pref.q))
    # println(relative_maximum_error(sol.p, pref.p))
    # println()
    @test relative_maximum_error(sol.q, pref.q) < 4E-8
    # @test relative_maximum_error(sol.p, pref.p) < 4E-8

    sol = integrate(lode, VPRKLobattoIIIBIIIB̄(4))
    # println(relative_maximum_error(sol.q, pref.q))
    # println(relative_maximum_error(sol.p, pref.p))
    # println()
    @test relative_maximum_error(sol.q, pref.q) < 8E-13
    # @test relative_maximum_error(sol.p, pref.p) < 8E-13

end
