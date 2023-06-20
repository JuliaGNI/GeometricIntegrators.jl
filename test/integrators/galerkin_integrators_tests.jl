using GeometricIntegrators
using GeometricEquations.Tests.HarmonicOscillator
using Test

using CompactBasisFunctions
using QuadratureRules

using GeometricEquations.Tests.HarmonicOscillator: reference_solution

iode = iodeproblem()

QGau4 = GaussLegendreQuadrature(4)
BGau4 = Lagrange(QuadratureRules.nodes(QGau4))


### CGVI Integrators ###

cgint = IntegratorCGVI(iode, BGau4, QGau4)
cgsol = integrate(iode, cgint)
@test relative_maximum_error(cgsol.q, reference_solution) < 1E-7


### DGVI Integrators ###

dgint = IntegratorDGVI(iode, BGau4, QGau4)
dgsol = integrate(iode, dgint)
@test relative_maximum_error(dgsol.q, reference_solution) < 1E-7

dgint = IntegratorDGVIP0(iode, BGau4, QGau4)
dgsol = integrate(iode, dgint)
@test relative_maximum_error(dgsol.q, reference_solution) < 1E-7

dgint = IntegratorDGVIP1(iode, BGau4, QGau4)
dgsol = integrate(iode, dgint)
@test relative_maximum_error(dgsol.q, reference_solution) < 1E-7

dgint = IntegratorDGVIEXP(iode, BGau4, QGau4)
dgsol = integrate(iode, dgint)
@test relative_maximum_error(dgsol.q, reference_solution) < 1E-7

dgint = IntegratorDGVIPI(iode, BGau4, QGau4, Discontinuity(PathIntegralLinear(), LobattoLegendreQuadrature(2)))
dgsol = integrate(iode, dgint)
@test relative_maximum_error(dgsol.q, reference_solution) < 1E-7
