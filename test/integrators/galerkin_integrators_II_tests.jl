using GeometricIntegrators
using GeometricProblems
using CompactBasisFunctions
using QuadratureRules
using Test


iode = GeometricProblems.HarmonicOscillator.iodeproblem()
pref = GeometricProblems.HarmonicOscillator.exact_solution(GeometricProblems.HarmonicOscillator.podeproblem())

QGau = GaussLegendreQuadrature(8)
BGau = Lagrange(QuadratureRules.nodes(QGau))


### CGVI Integrators ###
cgsol = integrate(iode, CGVI_II(BGau, QGau))
@test relative_maximum_error(cgsol.q, pref.q) < 1E-7

