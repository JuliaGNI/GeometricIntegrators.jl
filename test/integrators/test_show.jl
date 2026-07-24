using GeometricIntegrators
using GeometricProblems.HarmonicOscillator

ode  = odeproblem()
pode = podeproblem()
hode = hodeproblem()
iode = iodeproblem()
lode = lodeproblem()
sode = sodeproblem()
ldae = ldaeproblem()


show(stdout, GeometricIntegrator(ode, ExplicitEuler()))
show(stdout, GeometricIntegrator(ode, TableauCrouzeix()))
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

# FLRK (Formal Lagrangian RK) is not available in the current architecture
# (method commented out in method_list.jl; integrator source removed):
# show(stdout, GeometricIntegrator(lode, FLRK(Gauss(1))))

show(stdout, GeometricIntegrator(lode, DVIA()))
show(stdout, GeometricIntegrator(lode, DVIB()))
show(stdout, GeometricIntegrator(lode, CMDVI()))
show(stdout, GeometricIntegrator(lode, CTDVI()))

show(stdout, GeometricIntegrator(iode, DVRK(Gauss(1))))
