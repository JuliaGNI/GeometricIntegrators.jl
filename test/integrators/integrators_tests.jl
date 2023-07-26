using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
using Test


ode = odeproblem()
ref = exact_solution(ode)

@testset "$(rpad("Various integrators", 80))" begin

    sol = integrate(ode, ExplicitEuler())
    err = relative_maximum_error(sol, ref)
    @test err.q < 5E-2

    sol = integrate(ode, ImplicitEuler())
    err = relative_maximum_error(sol, ref)
    @test err.q < 5E-2

end
