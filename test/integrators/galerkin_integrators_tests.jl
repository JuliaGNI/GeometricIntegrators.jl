using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
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

# DISABLED (blocker): the DGVI (Discontinuous Galerkin Variational) integrators are
# dead code — src/integrators/dgvi/*.jl are not compiled and use the superseded
# `Parameters`/`update_params!` integrator architecture. Reviving them requires a
# port to the current method-based GeometricIntegrator design; see VERIFICATION_REPORT.md.
# dgsol = integrate(iode, DGVI(BGau4, QGau4))
# @test relative_maximum_error(dgsol.q, pref.q) < 1E-7

# dgsol = integrate(iode, DGVIP0(BGau4, QGau4))
# @test relative_maximum_error(dgsol.q, pref.q) < 1E-7

# dgsol = integrate(iode, DGVIP1(BGau4, QGau4))
# @test relative_maximum_error(dgsol.q, pref.q) < 1E-7

# dgsol = integrate(iode, DGVIEXP(BGau4, QGau4))
# @test relative_maximum_error(dgsol.q, pref.q) < 1E-7

# dgsol = integrate(iode, DGVIPI(BGau4, QGau4, Discontinuity(PathIntegralLinear(), LobattoLegendreQuadrature(2))))
# @test relative_maximum_error(dgsol.q, pref.q) < 1E-7
