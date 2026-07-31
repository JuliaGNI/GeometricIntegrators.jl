using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
import GeometricProblems.LotkaVolterra2d as LotkaVolterra2d
using CompactBasisFunctions
using QuadratureRules
using Test


iode = iodeproblem()
pref = exact_solution(podeproblem())

QGau4 = GaussLegendreQuadrature(4)
BGau4 = Lagrange(QuadratureRules.nodes(QGau4))


### CGVI Integrators ###

cgsol = integrate(iode, CGVI(BGau4, QGau4))
@test relative_maximum_error(cgsol.q, pref.q) < 8E-13


### DGVI Integrators ###

# DGVIs discretise a *fully degenerate* Lagrangian L = ϑ(q)·q̇ − H(q); `ϑ` is evaluated as
# `ϑ(θ, t, q, q)`, which is only meaningful when it does not depend on the velocity. On the
# regular harmonic oscillator used above for CGVI the closure equation of `DGVI` collapses
# to `q − p`, which is independent of every unknown and gives a singular Jacobian — which
# is why the previously-disabled tests here could never have run. The Lotka–Volterra
# `iodeproblem_dg` is provided by GeometricProblems for exactly this purpose.

const dg_q₀ = [1.0, 1.0]
const dg_params = (a₁=1.0, a₂=1.0, b₁=-1.0, b₂=-2.0)
const dg_tspan = (0.0, 1.0)
const dg_Δt = 0.05

dgiode = LotkaVolterra2d.iodeproblem_dg(dg_q₀; timespan=dg_tspan, timestep=dg_Δt, parameters=dg_params)
dgref = integrate(LotkaVolterra2d.odeproblem(dg_q₀; timespan=dg_tspan, timestep=dg_Δt, parameters=dg_params), Gauss(8))

dgjump = Discontinuity(PathIntegralLinear(), LobattoLegendreQuadrature(2))

@testset "$(rpad("Discontinuous Galerkin variational integrators",80))" begin

    # Tolerances measured at Δt = 0.05 against a Gauss(8) reference. The two variants with
    # an *average*-based flux (DGVIPI with a linear path, DGVIEXP with the midpoint
    # average) are some seven orders of magnitude more accurate than the three that
    # evaluate the flux at the nodal value — see
    # test/verification/dgvi_convergence_tests.jl for the order measurements behind that.
    #
    #            s = 4 measured
    #   DGVI     2.6E-8
    #   DGVIP0   2.6E-8
    #   DGVIP1   2.6E-8
    #   DGVIPI   8.0E-14
    #   DGVIEXP  5.2E-14

    dgsol = integrate(dgiode, DGVI(BGau4, QGau4))
    @test relative_maximum_error(dgsol.q, dgref.q) < 4E-8

    dgsol = integrate(dgiode, DGVIP0(BGau4, QGau4))
    @test relative_maximum_error(dgsol.q, dgref.q) < 4E-8

    dgsol = integrate(dgiode, DGVIP1(BGau4, QGau4))
    @test relative_maximum_error(dgsol.q, dgref.q) < 4E-8

    dgsol = integrate(dgiode, DGVIEXP(BGau4, QGau4))
    @test relative_maximum_error(dgsol.q, dgref.q) < 8E-14

    dgsol = integrate(dgiode, DGVIPI(BGau4, QGau4, dgjump))
    @test relative_maximum_error(dgsol.q, dgref.q) < 1E-13

    # `DGVIP1` adds a continuity constraint and a final projection on top of `DGVI`'s
    # equations; on a consistent solution the two agree to round-off.
    @test relative_maximum_error(integrate(dgiode, DGVI(BGau4, QGau4)).q,
                                 integrate(dgiode, DGVIP1(BGau4, QGau4)).q) < 1E-12

    # traits
    for m in (DGVI(BGau4, QGau4), DGVIP0(BGau4, QGau4), DGVIP1(BGau4, QGau4),
              DGVIEXP(BGau4, QGau4), DGVIPI(BGau4, QGau4, dgjump))
        @test isimplicit(m)
        @test !isexplicit(m)
        @test islodemethod(m)
        @test isiodemethod(m)
    end

    # the number of degrees of freedom the step adds on top of the basis coefficients
    @test GeometricIntegrators.Integrators.nclosure(DGVI(BGau4, QGau4)) == 1
    @test GeometricIntegrators.Integrators.nclosure(DGVIP0(BGau4, QGau4)) == 1
    @test GeometricIntegrators.Integrators.nclosure(DGVIPI(BGau4, QGau4, dgjump)) == 1
    @test GeometricIntegrators.Integrators.nclosure(DGVIP1(BGau4, QGau4)) == 2
    @test GeometricIntegrators.Integrators.nclosure(DGVIEXP(BGau4, QGau4)) == 2

    # `DGVIP0` contracts the basis' boundary coefficients with the one-form sampled at the
    # quadrature nodes, so it needs the two node sets to coincide. Enforced at
    # construction: without the check this combination builds fine and then raises a bare
    # `BoundsError` from inside the first step.
    @test_throws AssertionError DGVIP0(Lagrange(QuadratureRules.nodes(GaussLegendreQuadrature(3))), QGau4)

    # the other four variants place no such requirement on the pair
    let B3 = Lagrange(QuadratureRules.nodes(GaussLegendreQuadrature(3)))
        @test GeometricIntegrators.Integrators.nbasis(DGVI(B3, QGau4)) == 3
        @test GeometricIntegrators.Integrators.nnodes(DGVI(B3, QGau4)) == 4
    end

    # only the path-integral variant carries jump-quadrature coefficients
    @test GeometricIntegrators.Integrators.njump(DGVI(BGau4, QGau4)) == 0
    @test GeometricIntegrators.Integrators.njump(DGVIPI(BGau4, QGau4, dgjump)) == 2

    # `DGVIPI` with a linear path and 2-point Lobatto quadrature reduces to the
    # trapezoidal flux on the full jump, with the nodal value at the path midpoint
    @test DGVIPI(BGau4, QGau4, dgjump).ρ⁻ ≈ 0.5
    @test DGVIPI(BGau4, QGau4, dgjump).ρ⁺ ≈ 0.5

end
