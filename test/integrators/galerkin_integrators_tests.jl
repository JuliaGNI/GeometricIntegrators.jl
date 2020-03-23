
using GeometricIntegrators.BasisFunctions
using GeometricIntegrators.CommonFunctions
using GeometricIntegrators.Config
using GeometricIntegrators.Discontinuities
using GeometricIntegrators.Integrators
using GeometricIntegrators.Quadratures
using GeometricIntegrators.Solvers
using GeometricIntegrators.TestProblems.HarmonicOscillatorProblem
using GeometricIntegrators.Utils
using Test

using GeometricIntegrators.TestProblems.HarmonicOscillatorProblem: Δt, nt, refx, refq, refp

iode = harmonic_oscillator_iode()

QGau4 = GaussLegendreQuadrature(4)
BGau4 = LagrangeBasis(nodes(QGau4))


### CGVI Integrators ###

cgint = IntegratorCGVI(iode, BGau4, QGau4, Δt)
cgsol = integrate(iode, cgint, nt)
@test rel_err(cgsol.q, refx) < 1E-7


### DGVI Integrators ###

dgint = IntegratorDGVI(iode, BGau4, QGau4, Δt)
dgsol = integrate(iode, dgint, nt)
@test rel_err(dgsol.q, refx) < 1E-7

dgint = IntegratorDGVIP0(iode, BGau4, QGau4, Δt)
dgsol = integrate(iode, dgint, nt)
@test rel_err(dgsol.q, refx) < 1E-7

dgint = IntegratorDGVIP1(iode, BGau4, QGau4, Δt)
dgsol = integrate(iode, dgint, nt)
@test rel_err(dgsol.q, refx) < 1E-7

dgint = IntegratorDGVIEXP(iode, BGau4, QGau4, Δt)
dgsol = integrate(iode, dgint, nt)
@test rel_err(dgsol.q, refx) < 1E-7

dgint = IntegratorDGVIPI(iode, BGau4, QGau4, Discontinuity(PathIntegralLinear(), LobattoLegendreQuadrature(2)), Δt)
dgsol = integrate(iode, dgint, nt)
@test rel_err(dgsol.q, refx) < 1E-7
