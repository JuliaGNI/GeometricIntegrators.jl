
using GeometricBase
using GeometricIntegrators.Config
using GeometricIntegrators.Discontinuities
using GeometricIntegrators.Integrators
using GeometricIntegrators.Utils
using GeometricProblems.HarmonicOscillator
using Test

using CompactBasisFunctions
using QuadratureRules

using GeometricProblems.HarmonicOscillator: Δt, nt, refx, refq, refp

iode = harmonic_oscillator_iode()

QGau4 = GaussLegendreQuadrature(4)
BGau4 = Lagrange(nodes(QGau4))


### CGVI Integrators ###

cgint = IntegratorCGVI(iode, BGau4, QGau4, Δt)
cgsol = integrate(iode, cgint, nt)
@test relative_maximum_error(cgsol.q, refx) < 1E-7


### DGVI Integrators ###

dgint = IntegratorDGVI(iode, BGau4, QGau4, Δt)
dgsol = integrate(iode, dgint, nt)
@test relative_maximum_error(dgsol.q, refx) < 1E-7

dgint = IntegratorDGVIP0(iode, BGau4, QGau4, Δt)
dgsol = integrate(iode, dgint, nt)
@test relative_maximum_error(dgsol.q, refx) < 1E-7

dgint = IntegratorDGVIP1(iode, BGau4, QGau4, Δt)
dgsol = integrate(iode, dgint, nt)
@test relative_maximum_error(dgsol.q, refx) < 1E-7

dgint = IntegratorDGVIEXP(iode, BGau4, QGau4, Δt)
dgsol = integrate(iode, dgint, nt)
@test relative_maximum_error(dgsol.q, refx) < 1E-7

dgint = IntegratorDGVIPI(iode, BGau4, QGau4, Discontinuity(PathIntegralLinear(), LobattoLegendreQuadrature(2)), Δt)
dgsol = integrate(iode, dgint, nt)
@test relative_maximum_error(dgsol.q, refx) < 1E-7
