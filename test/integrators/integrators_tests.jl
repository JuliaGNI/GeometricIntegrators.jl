using GeometricIntegrators
using SimpleSolvers
using Test

using GeometricProblems.HarmonicOscillator
using GeometricProblems.HarmonicOscillator: reference_solution

ode  = odeproblem()


@testset "$(rpad("Various integrators",80))" begin

    sol = integrate(ode, ExplicitEuler())
    @test relative_maximum_error(sol.q, reference_solution) < 5E-2

    sol = integrate(ode, ImplicitEuler())
    @test relative_maximum_error(sol.q, reference_solution) < 5E-2

end
