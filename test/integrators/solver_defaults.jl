using GeometricIntegrators
using GeometricIntegrators.SPARK
using InteractiveUtils: subtypes
using SimpleSolvers
using Test

import GeometricIntegratorsBase
import GeometricIntegratorsBase: solver
import GeometricProblems.HarmonicOscillator
import GeometricProblems.LotkaVolterra2d

const GI = GeometricIntegrators

linearmethod(s) = SimpleSolvers.method(SimpleSolvers.linearsolver(s))

ode = HarmonicOscillator.odeproblem()
dae = HarmonicOscillator.daeproblem()
hdae = LotkaVolterra2d.hdaeproblem()

@testset "$(rpad("Linear solver of the nonlinear solvers",80))" begin
    # the GeometricIntegratorsBase `initsolver` path, for Newton, QuasiNewton and DogLeg
    for s in (Newton(), QuasiNewton(), DogLeg())
        @test linearmethod(solver(GeometricIntegrator(ode, Gauss(2); solver = s))) isa
              SimpleSolvers.LU
    end
    @test linearmethod(solver(GeometricIntegrator(hdae, TableauHSPARKLobattoIIIAB(2)))) isa
          SimpleSolvers.LU

    # the caller's value wins
    int = GeometricIntegrator(ode, Gauss(2); linear_solver_method = SimpleSolvers.LapackLU())
    @test linearmethod(solver(int)) isa SimpleSolvers.LapackLU

    # the DIRK stage solvers
    int = GeometricIntegrator(ode, Crouzeix())
    @test all(linearmethod(s) isa SimpleSolvers.LU for s in solver(int).solvers)
    int = GeometricIntegrator(ode, Crouzeix(); linear_solver_method = SimpleSolvers.LapackLU())
    @test all(linearmethod(s) isa SimpleSolvers.LapackLU for s in solver(int).solvers)

    # the projection solver and the solver of the projected method
    int = GeometricIntegrator(dae, PostProjection(Gauss(1)))
    @test linearmethod(solver(int)) isa SimpleSolvers.LU
    @test linearmethod(solver(int.subint)) isa SimpleSolvers.LU
    int = GeometricIntegrator(dae, MidpointProjection(Gauss(1)))
    @test linearmethod(solver(int)) isa SimpleSolvers.LU

    # options without `linear_solver_method` replace the projection defaults, not the LU default
    int = GeometricIntegrator(dae, PostProjection(Gauss(1)); f_abstol = 1e-14)
    @test linearmethod(solver(int)) isa SimpleSolvers.LU
    int = GeometricIntegrator(dae, MidpointProjection(Gauss(1)); f_abstol = 1e-14)
    @test linearmethod(solver(int)) isa SimpleSolvers.LU
    int = GeometricIntegrator(dae, PostProjection(Gauss(1));
        linear_solver_method = SimpleSolvers.LapackLU())
    @test linearmethod(solver(int)) isa SimpleSolvers.LapackLU

    # Picard takes no `linear_solver_method`
    @test solver(GeometricIntegrator(ode, Gauss(2); solver = Picard())) isa
          SimpleSolvers.NonlinearSolver
end

@testset "$(rpad("Every method of this package gets the LU default",80))" begin
    own(T) = Base.moduleroot(parentmodule(T)) === GI
    roots = Any[]
    walk(T) = foreach(S -> own(S) ? push!(roots, S) : walk(S), subtypes(T))
    walk(GeometricIntegratorsBase.GeometricMethod)
    @test !isempty(roots)
    for T in roots, S in (Newton, DogLeg)

        m = which(GeometricIntegratorsBase.initsolver,
            Tuple{S, T, GeometricIntegratorsBase.CacheDict})
        @test Base.moduleroot(m.module) === GI
    end
end
