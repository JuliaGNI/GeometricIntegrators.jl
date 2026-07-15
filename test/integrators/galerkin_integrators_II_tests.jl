using GeometricIntegrators
using GeometricProblems
using CompactBasisFunctions
using QuadratureRules
using Test


iode = GeometricProblems.HarmonicOscillator.iodeproblem()
pref = GeometricProblems.HarmonicOscillator.exact_solution(GeometricProblems.HarmonicOscillator.podeproblem())

# Use the Lobatto Legendre Quadrature Rule only!
QGau = LobattoLegendreQuadrature(4)
BGau = Lagrange(QuadratureRules.nodes(QGau))

### CGVI Integrators ###
cgsol_ii = integrate(iode, CGVI_II(BGau, QGau))
@test relative_maximum_error(cgsol_ii.q, pref.q) < 1E-7



