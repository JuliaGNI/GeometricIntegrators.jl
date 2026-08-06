using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricProblems.HarmonicOscillator
import GeometricProblems.LotkaVolterra2d as LotkaVolterra2d
import GeometricProblems.LotkaVolterra2dSingular
using CompactBasisFunctions
using QuadratureRules
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

# The Galerkin methods need a basis and a quadrature rule, the DGVIs additionally a
# fully degenerate Lagrangian, and the SPARK methods a DAE formulation — the same
# constructions as in test/integrators/galerkin_integrators_tests.jl and
# test/spark/spark_integrators_tests.jl.
QGau4 = GaussLegendreQuadrature(4)
BGau4 = Lagrange(QuadratureRules.nodes(QGau4))
dgjump = Discontinuity(PathIntegralLinear(), LobattoLegendreQuadrature(2))

const lv_q₀ = [1.0, 1.0]
const lv_params = (a₁=1.0, a₂=1.0, b₁=-1.0, b₂=-2.0)
const lv_tspan = (0.0, 0.1)
const lv_Δt = 0.01

lvargs = (timespan=lv_tspan, timestep=lv_Δt, parameters=lv_params)

dgiode = LotkaVolterra2d.iodeproblem_dg(lv_q₀; lvargs...)
lvidae = LotkaVolterra2d.idaeproblem(lv_q₀; lvargs...)
lvldae = LotkaVolterra2d.ldaeproblem(lv_q₀; lvargs...)
lvpdae = LotkaVolterra2d.pdaeproblem(lv_q₀; lvargs...)
lvhdae = LotkaVolterra2d.hdaeproblem(lv_q₀; lvargs...)
lvslrk = LotkaVolterra2d.ldaeproblem_slrk(lv_q₀; lvargs...)

# The `show` output is captured with `sprint` rather than written to `stdout`: it is
# verbose (tableaus, coefficients, literature references) and would otherwise bury the
# test summary. Where the object carries a `show` method of its own, its title line is
# asserted, so that a method printing nothing — or Julia's default struct dump silently
# standing in for a missing one — fails instead of passing.
test_show(x, title) = @test occursin(title, sprint(show, x))

# For objects without a specialised `show` all that can be checked is that the default
# dump runs and produces output. They are listed separately below so the gaps stay
# visible rather than hiding behind a non-emptiness check.
test_show_default(x) = @test length(sprint(show, x)) > 0


@testset "$(rpad("Show Methods and Integrators",80))" begin

    ### the methods themselves ###

    # `RKMethod` and `PRKMethod` (src/integrators/rk/abstract.jl), `VPRKMethod`
    # (src/integrators/vi/vprk_methods.jl), `CGVI` (src/integrators/cgvi/integrators_cgvi.jl)
    # and `DGVIMethod` (src/integrators/dgvi/integrators_dgvi_common.jl) print their tableau
    # or basis; the remaining method structs have no `show` of their own.
    test_show(RK(TableauCrouzeix()), "Runge-Kutta Method with Tableau")
    test_show(Gauss(1), "Runge-Kutta Method with Tableau")
    test_show(LobattoIIIAIIIB(2), "Partitioned Runge-Kutta Method with Tableau")
    test_show(VPRK(Gauss(1)), "Variational Partitioned Runge-Kutta Method with Tableau")
    test_show(CGVI(BGau4, QGau4), "Continuous Galerkin Variational Integrator")
    test_show(DGVI(BGau4, QGau4), "Discontinuous Galerkin Variational Integrator")

    test_show_default(ExplicitEuler())
    test_show_default(VPRKpMidpoint(Gauss(1)))


    ### the integrators ###

    test_show(GeometricIntegrator(ode, RK(TableauCrouzeix())), "Diagonally Implicit Runge-Kutta Integrator")
    test_show(GeometricIntegrator(ode, Gauss(1)), "Implicit Runge-Kutta Integrator")
    test_show(GeometricIntegrator(ode, PGLRK(3)), "Projected Gauss-Legendre Runge-Kutta Integrator")

    test_show(GeometricIntegrator(pode, LobattoIIIAIIIB(2)), "Explicit Partitioned Runge-Kutta Integrator")
    test_show(GeometricIntegrator(pode, Gauss(1)), "Implicit Partitioned Runge-Kutta Integrator")

    test_show(GeometricIntegrator(hode, SymplecticEulerA()), "Explicit Partitioned Runge-Kutta Integrator")
    test_show(GeometricIntegrator(hode, Gauss(1)), "Implicit Partitioned Runge-Kutta Integrator")

    test_show(GeometricIntegrator(iode, Gauss(1)), "Runge-Kutta Integrator for Implicit Equations")
    test_show(GeometricIntegrator(iode, LobattoIIIAIIIB(2)), "Explicit Partitioned Runge-Kutta Integrator")
    test_show(GeometricIntegrator(iode, VPRKpTableau(4)), "Variational Partitioned Runge-Kutta Integrator with Projection in the Tableau")

    test_show(GeometricIntegrator(lode, FLRK(Gauss(1))), "Formal Lagrangian Runge-Kutta Integrator")

    test_show(GeometricIntegrator(lode, PMVImidpoint()), "Midpoint variational integrator in position-momentum form")
    test_show(GeometricIntegrator(lode, PMVItrapezoidal()), "Trapezoidal variational integrator in position-momentum form")

    # the Hamilton-Pontryagin integrators take the discrete velocity map and its
    # derivatives, as in test/integrators/hamilton_pontryagin_integrators_tests.jl
    ϕ(v, q̄, q, a, Δt) = v .= (q .- q̄) ./ Δt

    function D₁ϕ(d, q̄, q, a, Δt)
        d .= 0
        for i in eachindex(q)
            d[i, i] = -1 / Δt
        end
    end

    function D₂ϕ(d, q̄, q, a, Δt)
        d .= 0
        for i in eachindex(q)
            d[i, i] = +1 / Δt
        end
    end

    Dₐϕ(d, q̄, q, a, Δt) = d .= 0

    test_show(GeometricIntegrator(lode, HPImidpoint(ϕ, D₁ϕ, D₂ϕ, Dₐϕ, Float64[])), "Hamilton-Pontryagin Integrator using midpoint quadrature")
    test_show(GeometricIntegrator(lode, HPItrapezoidal(ϕ, D₁ϕ, D₂ϕ, Dₐϕ, Float64[])), "Hamilton-Pontryagin Integrator using trapezoidal quadrature")

    test_show(GeometricIntegrator(dgiode, DGVI(BGau4, QGau4)), "Discontinuous Galerkin Variational Integrator")
    test_show(GeometricIntegrator(dgiode, DGVIPI(BGau4, QGau4, dgjump)), "Discontinuous Galerkin Variational Integrator")

    test_show(GeometricIntegrator(dvilode, DVIA()), "Degenerate Variational Integrator (Euler-A)")
    test_show(GeometricIntegrator(dvilode, DVIB()), "Degenerate Variational Integrator (Euler-B)")
    test_show(GeometricIntegrator(dvilode, CMDVI()), "Midpoint Degenerate Variational Integrator")
    test_show(GeometricIntegrator(dvilode, CTDVI()), "Trapezoidal Degenerate Variational Integrator")
    test_show(GeometricIntegrator(dvilode, DVRK(Gauss(1))), "Runge-Kutta Integrator for Degenerate Lagrangians")


    ### the SPARK integrators ###

    # only constructed, never integrated, so none of the solver warnings the SPARK tests
    # silence can appear here. All ten print the tableau via `tableau(method(int)).q/.p`;
    # spelling that `method(int).q` threw a `FieldError` for the seven method types that
    # wrap their tableau instead of being one.
    test_show(GeometricIntegrator(lvslrk, SLRKLobattoIIIAB(2)), "Specialised Partitioned Additive Runge-Kutta integrator for degenerate")
    test_show(GeometricIntegrator(lvidae, SPARKGLRK(1)), "Specialised Partitioned Additive Runge-Kutta integrator for index-two DAE systems")
    test_show(GeometricIntegrator(lvidae, TableauGausspSymplectic(1)), "Variational partitioned additive Runge-Kutta integrator")
    test_show(GeometricIntegrator(lvidae, VSPARK(SPARKLobABC(3))), "Specialised Partitioned Additive Runge-Kutta integrator for Variational systems")
    test_show(GeometricIntegrator(lvidae, TableauVSPARKGLRKpMidpoint(1)), "with projection on primary constraint")
    test_show(GeometricIntegrator(lvldae, TableauVSPARKLobattoIIIAB(2)), "variational systems with projection on secondary constraint")
    test_show(GeometricIntegrator(lvpdae, TableauHPARKGLRK(1)), "Partitioned Additive Runge-Kutta integrator for Hamiltonian systems subject")
    test_show(GeometricIntegrator(lvpdae, HSPARK(SPARKGLRK(1))), "Specialised Partitioned Additive Runge-Kutta integrator for Hamiltonian systems *EXPERIMENTAL*")
    test_show(GeometricIntegrator(lvpdae, TableauHSPARKGLRKpSymmetric(1)), "with projection on primary constraint *EXPERIMENTAL*")
    test_show(GeometricIntegrator(lvhdae, TableauHSPARKLobattoIIIAB(2)), "subject to Dirac constraints with projection on secondary constraint")


    ### integrators that fall back to Julia's default struct dump ###

    # `ExplicitEuler` and `ImplicitEuler` are implemented in GeometricIntegratorsBase, and
    # the `{<:ERK}` / `{<:IRK}` methods in src/integrators/rk/ do not cover them, because
    # neither method is converted to `ERK`/`IRK` the way `Gauss` is.
    test_show_default(GeometricIntegrator(ode, ExplicitEuler()))
    test_show_default(GeometricIntegrator(ode, ImplicitEuler()))

    # the splitting and composition integrators define no `show` at all
    test_show_default(GeometricIntegrator(sode, LieA()))
    test_show_default(GeometricIntegrator(sode, Composition(Strang())))

    # `show` commented out in src/integrators/vprk/integrators_vprk_abstract.jl
    test_show_default(GeometricIntegrator(iode, VPRK(Gauss(1))))

    # the projected integrators: `show` commented out in
    # src/projections/{midpoint,standard,symmetric}_projection.jl
    test_show_default(GeometricIntegrator(iode, VPRKpInternal(Gauss(1))))
    test_show_default(GeometricIntegrator(iode, VPRKpMidpoint(Gauss(1))))
    test_show_default(GeometricIntegrator(ldae, VPRKpSecondary(Gauss(1))))
    test_show_default(GeometricIntegrator(iode, VPRKpStandard(Gauss(1))))
    test_show_default(GeometricIntegrator(iode, VPRKpSymmetric(Gauss(1))))
    test_show_default(GeometricIntegrator(iode, VPRKpVariational(Gauss(1))))

    # `CGVI` defines `show` for the method (asserted above), not for the integrator
    test_show_default(GeometricIntegrator(iode, CGVI(BGau4, QGau4)))

end
