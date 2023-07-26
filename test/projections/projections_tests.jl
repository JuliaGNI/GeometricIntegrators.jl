using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
using Test


dae = daeproblem()
ref = exact_solution(odeproblem())


@testset "$(rpad("Post-projection with Runge-Kutta integrators",80))" begin

    sol = integrate(dae, PostProjection(RK4()))
    @test relative_maximum_error(sol.q, ref.q) < 2E-7
    @test abs(hamiltonian(sol.t[end], sol.q[end], dae.parameters) - hamiltonian(sol.t[begin], sol.q[begin], dae.parameters)) < eps()

    sol = integrate(dae, PostProjection(Gauss(1)))
    @test relative_maximum_error(sol.q, ref.q) < 4E-4
    @test abs(hamiltonian(sol.t[end], sol.q[end], dae.parameters) - hamiltonian(sol.t[begin], sol.q[begin], dae.parameters)) < eps()

    sol = integrate(dae, PostProjection(Gauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 4E-8
    @test abs(hamiltonian(sol.t[end], sol.q[end], dae.parameters) - hamiltonian(sol.t[begin], sol.q[begin], dae.parameters)) < eps()

    sol = integrate(dae, PostProjection(Gauss(3)))
    @test relative_maximum_error(sol.q, ref.q) < 2E-12
    @test abs(hamiltonian(sol.t[end], sol.q[end], dae.parameters) - hamiltonian(sol.t[begin], sol.q[begin], dae.parameters))  < eps()

    sol = integrate(dae, PostProjection(Gauss(4)))
    @test relative_maximum_error(sol.q, ref.q) < 4E-15
    @test abs(hamiltonian(sol.t[end], sol.q[end], dae.parameters) - hamiltonian(sol.t[begin], sol.q[begin], dae.parameters)) < eps()

end


@testset "$(rpad("Midpoint projection with Runge-Kutta integrators",80))" begin

    # sol = integrate(dae, MidpointProjection(RK4()))
    # @test relative_maximum_error(sol.q, ref.q) < 2E-7
    # @test abs(hamiltonian(sol.t[end], sol.q[end], dae.parameters) - hamiltonian(sol.t[begin], sol.q[begin], dae.parameters)) < eps()

    sol = integrate(dae, MidpointProjection(Gauss(1)))
    @test relative_maximum_error(sol.q, ref.q) < 4E-4
    @test abs(hamiltonian(sol.t[end], sol.q[end], dae.parameters) - hamiltonian(sol.t[begin], sol.q[begin], dae.parameters)) < 2eps()

    sol = integrate(dae, MidpointProjection(Gauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 4E-8
    @test abs(hamiltonian(sol.t[end], sol.q[end], dae.parameters) - hamiltonian(sol.t[begin], sol.q[begin], dae.parameters)) < 2eps()

    sol = integrate(dae, MidpointProjection(Gauss(3)))
    @test relative_maximum_error(sol.q, ref.q) < 2E-12
    @test abs(hamiltonian(sol.t[end], sol.q[end], dae.parameters) - hamiltonian(sol.t[begin], sol.q[begin], dae.parameters)) < 2eps()

    sol = integrate(dae, MidpointProjection(Gauss(4)))
    @test relative_maximum_error(sol.q, ref.q) < 4E-15
    @test abs(hamiltonian(sol.t[end], sol.q[end], dae.parameters) - hamiltonian(sol.t[begin], sol.q[begin], dae.parameters)) < 2eps()

end


@testset "$(rpad("Symmetric projection with Runge-Kutta integrators",80))" begin

    # sol = integrate(dae, SymmetricProjection(RK4()))
    # @test relative_maximum_error(sol.q, ref.q) < 2E-7
    # @test abs(hamiltonian(sol.t[end], sol.q[end], dae.parameters) - hamiltonian(sol.t[begin], sol.q[begin], dae.parameters)) < eps()

    sol = integrate(dae, SymmetricProjection(Gauss(1)))
    @test relative_maximum_error(sol.q, ref.q) < 4E-4
    @test abs(hamiltonian(sol.t[end], sol.q[end], dae.parameters) - hamiltonian(sol.t[begin], sol.q[begin], dae.parameters)) < 2eps()

    sol = integrate(dae, SymmetricProjection(Gauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 4E-8
    @test abs(hamiltonian(sol.t[end], sol.q[end], dae.parameters) - hamiltonian(sol.t[begin], sol.q[begin], dae.parameters)) < 2eps()

    sol = integrate(dae, SymmetricProjection(Gauss(3)))
    @test relative_maximum_error(sol.q, ref.q) < 2E-12
    @test abs(hamiltonian(sol.t[end], sol.q[end], dae.parameters) - hamiltonian(sol.t[begin], sol.q[begin], dae.parameters)) < 2eps()

    sol = integrate(dae, SymmetricProjection(Gauss(4)))
    @test relative_maximum_error(sol.q, ref.q) < 4E-15
    @test abs(hamiltonian(sol.t[end], sol.q[end], dae.parameters) - hamiltonian(sol.t[begin], sol.q[begin], dae.parameters)) < 2eps()

end
