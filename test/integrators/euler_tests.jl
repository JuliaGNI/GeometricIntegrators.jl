using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
using SimpleSolvers: Options
using Test


ode = odeproblem()
ref = exact_solution(ode)

@testset "$(rpad("Euler integrators", 80))" begin

    sol = integrate(ode, ExplicitEuler())
    err = relative_maximum_error(sol, ref)
    @test err.q < 5E-2

    sol = integrate(ode, ImplicitEuler())
    err = relative_maximum_error(sol, ref)
    @test err.q < 5E-2

    sol = integrate(ode, ImplicitEuler(); options = Options(
        min_iterations = 1,
        x_abstol = 2eps(),
        f_abstol = 2eps(),
    ))
    err = relative_maximum_error(sol, ref)
    @test err.q < 5E-2

end
