using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
import GeometricProblems.LotkaVolterra2dSingular
using RungeKutta.Tableaus: TableauCrouzeix
using Test

ode  = odeproblem()
pode = podeproblem()
hode = hodeproblem()
iode = iodeproblem()
lode = lodeproblem()
sode = sodeproblem()
ldae = ldaeproblem()

# The degenerate variational integrators need an even-dimensional configuration
# space: they split the coordinates into the half that carries the quadrature
# update and the half determined by the constraint p = ϑ(q). The harmonic
# oscillator problems used above have D = 1, for which none of these methods is
# defined and the cache constructors now throw, so use an in-class 2d problem here.
dvilode = LotkaVolterra2dSingular.lodeproblem()


# The `show` methods are captured with `sprint` rather than written to `stdout`:
# their output is verbose (tableaus, coefficients, literature references) and would
# otherwise bury the test summary. Asserting a non-empty result keeps the smoke test
# honest — a `show` method printing nothing at all now fails instead of passing
# silently.
@testset "$(rpad("Show Methods and Integrators",80))" begin

    # the methods themselves: one per family with a `show` method of its own
    # (`RKMethod` and `PRKMethod` in src/integrators/rk/abstract.jl, `VPRKMethod` in
    # src/integrators/vi/vprk_methods.jl). The Galerkin methods `CGVI` and `DGVIMethod`
    # also define `show`, but constructing them needs a basis and a quadrature rule;
    # they are left to the Galerkin tests.
    @test length(sprint(show, ExplicitEuler())) > 0
    @test length(sprint(show, RK(TableauCrouzeix()))) > 0
    @test length(sprint(show, Gauss(1))) > 0
    @test length(sprint(show, LobattoIIIAIIIB(2))) > 0
    @test length(sprint(show, VPRK(Gauss(1)))) > 0
    @test length(sprint(show, VPRKpMidpoint(Gauss(1)))) > 0

    # the integrators
    @test length(sprint(show, GeometricIntegrator(ode, ExplicitEuler()))) > 0
    @test length(sprint(show, GeometricIntegrator(ode, RK(TableauCrouzeix())))) > 0
    @test length(sprint(show, GeometricIntegrator(ode, ImplicitEuler()))) > 0

    @test length(sprint(show, GeometricIntegrator(pode, LobattoIIIAIIIB(2)))) > 0
    @test length(sprint(show, GeometricIntegrator(pode, Gauss(1)))) > 0

    @test length(sprint(show, GeometricIntegrator(hode, SymplecticEulerA()))) > 0
    @test length(sprint(show, GeometricIntegrator(hode, Gauss(1)))) > 0

    @test length(sprint(show, GeometricIntegrator(sode, LieA()))) > 0
    @test length(sprint(show, GeometricIntegrator(sode, Composition(Strang())))) > 0

    @test length(sprint(show, GeometricIntegrator(iode, Gauss(1)))) > 0
    @test length(sprint(show, GeometricIntegrator(iode, LobattoIIIAIIIB(2)))) > 0

    @test length(sprint(show, GeometricIntegrator(iode, VPRK(Gauss(1))))) > 0
    @test length(sprint(show, GeometricIntegrator(iode, VPRKpInternal(Gauss(1))))) > 0
    @test length(sprint(show, GeometricIntegrator(iode, VPRKpMidpoint(Gauss(1))))) > 0
    @test length(sprint(show, GeometricIntegrator(ldae, VPRKpSecondary(Gauss(1))))) > 0
    @test length(sprint(show, GeometricIntegrator(iode, VPRKpStandard(Gauss(1))))) > 0
    @test length(sprint(show, GeometricIntegrator(iode, VPRKpSymmetric(Gauss(1))))) > 0
    @test length(sprint(show, GeometricIntegrator(iode, VPRKpVariational(Gauss(1))))) > 0

    @test length(sprint(show, GeometricIntegrator(lode, FLRK(Gauss(1))))) > 0

    @test length(sprint(show, GeometricIntegrator(dvilode, DVIA()))) > 0
    @test length(sprint(show, GeometricIntegrator(dvilode, DVIB()))) > 0
    @test length(sprint(show, GeometricIntegrator(dvilode, CMDVI()))) > 0
    @test length(sprint(show, GeometricIntegrator(dvilode, CTDVI()))) > 0

    @test length(sprint(show, GeometricIntegrator(dvilode, DVRK(Gauss(1))))) > 0

end
