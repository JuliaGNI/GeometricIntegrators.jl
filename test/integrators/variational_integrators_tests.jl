using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
using Test

using GeometricEquations: DELEProblem


lode = lodeproblem()
pref = exact_solution(podeproblem())


@testset "$(rpad("Discrete Euler-Lagrange integrators",80))" begin

    sol = integrate(deleproblem_midpoint(), DiscreteEulerLagrange())
    @test relative_maximum_error(sol.q, pref.q) < 4E-4

    dele = DELEProblem(lode, Midpoint())
    dsol = integrate(dele, DiscreteEulerLagrange())
    @test relative_maximum_error(dsol.q, pref.q) < 4E-4
    @test relative_maximum_error(dsol.q, sol.q) < 1E-14

    sol = integrate(deleproblem_trapezoidal(), DiscreteEulerLagrange())
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
    @test relative_maximum_error(sol.q, pref.q) < 4E-4
    @test relative_maximum_error(sol.p, pref.p) < 8E-4

    ref = integrate(lode, PMVImidpoint())
    @test relative_maximum_error(sol.q, ref.q) < 8 * eps()

    sol = integrate(lode, VPRKGauss(2))
    @test relative_maximum_error(sol.q, pref.q) < 4E-8
    @test relative_maximum_error(sol.p, pref.p) < 8E-8

    sol = integrate(lode, VPRKGauss(3))
    @test relative_maximum_error(sol.q, pref.q) < 8E-13
    @test relative_maximum_error(sol.p, pref.p) < 4E-12

    sol = integrate(lode, VPRKGauss(4))
    @test relative_maximum_error(sol.q, pref.q) < 4E-16
    @test relative_maximum_error(sol.p, pref.p) < 1E-15


    sol = integrate(lode, VPRKLobattoIIIAIIIĀ(2))
    @test relative_maximum_error(sol.q, pref.q) < 2E-4
    @test relative_maximum_error(sol.p, pref.p) < 1E-3

    ref = integrate(lode, PMVItrapezoidal())
    @test relative_maximum_error(sol.q, ref.q) < 8 * eps()

    sol = integrate(lode, VPRKLobattoIIIAIIIĀ(3))
    @test relative_maximum_error(sol.q, pref.q) < 8E-9
    @test relative_maximum_error(sol.p, pref.p) < 2E-7

    sol = integrate(lode, VPRKLobattoIIIAIIIĀ(4))
    @test relative_maximum_error(sol.q, pref.q) < 2E-13
    @test relative_maximum_error(sol.p, pref.p) < 4E-12


    sol = integrate(lode, VPRKLobattoIIIBIIIB̄(2))
    @test relative_maximum_error(sol.q, pref.q) < 4E-4
    @test relative_maximum_error(sol.p, pref.p) < 8E-4

    sol = integrate(lode, VPRKLobattoIIIBIIIB̄(3))
    @test relative_maximum_error(sol.q, pref.q) < 4E-8
    @test relative_maximum_error(sol.p, pref.p) < 8E-8

    sol = integrate(lode, VPRKLobattoIIIBIIIB̄(4))
    @test relative_maximum_error(sol.q, pref.q) < 8E-13
    @test relative_maximum_error(sol.p, pref.p) < 4E-12

end
