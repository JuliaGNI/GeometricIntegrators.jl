
using GeometricIntegrators.BasisFunctions
using GeometricIntegrators.CommonFunctions
using GeometricIntegrators.Config
using GeometricIntegrators.Integrators
using GeometricIntegrators.Quadratures
using GeometricIntegrators.Solvers
using GeometricIntegrators.TestProblems.HarmonicOscillator
using GeometricIntegrators.Utils
using Test

using GeometricIntegrators.TestProblems.HarmonicOscillator: Δt, nt, refx, refq, refp

iode = oscillator_iode()

QGau4 = GaussLegendreQuadrature(4)
BGau4 = LagrangeBasis(nodes(QGau4))


### CGVI Integrators ###

cgint = IntegratorCGVI(iode, BGau4, QGau4, Δt)
cgsol = integrate(cgint, nt)

@test rel_err(cgsol.q, refx) < 1E-7


### DGVI Integrators ###

dgint = IntegratorDGVI(iode, BGau4, QGau4, Δt)
dgsol = integrate(dgint, nt)

@test rel_err(dgsol.q, refx) < 1E-7
