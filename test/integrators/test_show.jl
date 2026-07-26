using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
import GeometricProblems.LotkaVolterra2dSingular
using RungeKutta.Tableaus: TableauCrouzeix

ode  = odeproblem()
pode = podeproblem()
hode = hodeproblem()
iode = iodeproblem()
lode = lodeproblem()
sode = sodeproblem()
ldae = ldaeproblem()


show(stdout, GeometricIntegrator(ode, ExplicitEuler()))
show(stdout, GeometricIntegrator(ode, RK(TableauCrouzeix())))
show(stdout, GeometricIntegrator(ode, ImplicitEuler()))

show(stdout, GeometricIntegrator(pode, LobattoIIIAIIIB(2)))
show(stdout, GeometricIntegrator(pode, Gauss(1)))

show(stdout, GeometricIntegrator(iode, Gauss(1)))
show(stdout, GeometricIntegrator(iode, LobattoIIIAIIIB(2)))

show(stdout, GeometricIntegrator(iode, VPRK(Gauss(1))))
show(stdout, GeometricIntegrator(iode, VPRKpInternal(Gauss(1))))
show(stdout, GeometricIntegrator(iode, VPRKpMidpoint(Gauss(1))))
show(stdout, GeometricIntegrator(ldae, VPRKpSecondary(Gauss(1))))
show(stdout, GeometricIntegrator(iode, VPRKpStandard(Gauss(1))))
show(stdout, GeometricIntegrator(iode, VPRKpSymmetric(Gauss(1))))
show(stdout, GeometricIntegrator(iode, VPRKpVariational(Gauss(1))))

show(stdout, GeometricIntegrator(lode, FLRK(Gauss(1))))

# The degenerate variational integrators need an even-dimensional configuration
# space: they split the coordinates into the half that carries the quadrature
# update and the half determined by the constraint p = ϑ(q). The harmonic
# oscillator problems used above have D = 1, for which none of these methods is
# defined and the cache constructors now throw, so use an in-class 2d problem here.
dvilode = LotkaVolterra2dSingular.lodeproblem()

show(stdout, GeometricIntegrator(dvilode, DVIA()))
show(stdout, GeometricIntegrator(dvilode, DVIB()))
show(stdout, GeometricIntegrator(dvilode, CMDVI()))
show(stdout, GeometricIntegrator(dvilode, CTDVI()))

show(stdout, GeometricIntegrator(dvilode, DVRK(Gauss(1))))
