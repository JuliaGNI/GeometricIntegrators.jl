using GeometricIntegrators
using Test

using GeometricProblems.HarmonicOscillator
using GeometricProblems.HarmonicOscillator: reference_solution, reference_solution_q, reference_solution_p

iode = iodeproblem()


@testset "$(rpad("Vartiational Partitioned Runge-Kutta integrators",80))" begin

    sol = integrate(iode, VPRKGauss(1))
    @test relative_maximum_error(sol.q, reference_solution_q) < 4E-4

    sol = integrate(iode, VPRKGauss(2))
    @test relative_maximum_error(sol.q, reference_solution_q) < 4E-8

    sol = integrate(iode, VPRKGauss(3))
    @test relative_maximum_error(sol.q, reference_solution_q) < 8E-13

    sol = integrate(iode, VPRKLobattoIIIAIIIĀ(2))
    @test relative_maximum_error(sol.q, reference_solution_q) < 2E-4

    sol = integrate(iode, VPRKLobattoIIIAIIIĀ(3))
    @test relative_maximum_error(sol.q, reference_solution_q) < 8E-9

    sol = integrate(iode, VPRKLobattoIIIAIIIĀ(4))
    @test relative_maximum_error(sol.q, reference_solution_q) < 2E-13

    sol = integrate(iode, VPRKLobattoIIIBIIIB̄(2))
    @test relative_maximum_error(sol.q, reference_solution_q) < 4E-4

    sol = integrate(iode, VPRKLobattoIIIBIIIB̄(3))
    @test relative_maximum_error(sol.q, reference_solution_q) < 4E-8

    sol = integrate(iode, VPRKLobattoIIIBIIIB̄(4))
    @test relative_maximum_error(sol.q, reference_solution_q) < 8E-13

end
