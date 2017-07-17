
set_config(:nls_solver, NewtonSolver)

Δt = 0.1
nt = 10

ode  = pendulum_ode()
pode = pendulum_pode()
iode = pendulum_iode()
idae = pendulum_idae()

glrk1 = getTableauGLRK(1)
glrk2 = getTableauGLRK(2)

@test typeof(Integrator(ode, getTableauExplicitMidpoint(), Δt)) <: IntegratorERK
@test typeof(Integrator(ode, getTableauCrouzeix(), Δt)) <: IntegratorDIRK
@test typeof(Integrator(ode, getTableauImplicitMidpoint(), Δt)) <: IntegratorFIRK

int = Integrator(ode, getTableauExplicitEuler(), Δt)
sol = integrate(int, nt)

int = Integrator(ode, getTableauImplicitEuler(), Δt)
sol = integrate(int, nt)

int = IntegratorFIRK(ode, getTableauImplicitMidpoint(), Δt)
sol = integrate(int, nt)

int = IntegratorFIRK(ode, getTableauGLRK(1), Δt)
sol = integrate(int, nt)

int = IntegratorFIRK(ode, getTableauGLRK(2), Δt)
sol = integrate(int, nt)

int = IntegratorFIRK(ode, getTableauGLRK(3), Δt)
sol = integrate(int, nt)

int = IntegratorFIRK(ode, getTableauGLRK(4), Δt)
sol = integrate(int, nt)

int = IntegratorFIRK(ode, getTableauGLRK(5), Δt)
sol = integrate(int, nt)

int = IntegratorFIRK(ode, getTableauGLRK(6), Δt)
sol = integrate(int, nt)

int = IntegratorFIRK(ode, getTableauGLRK(7), Δt)
sol = integrate(int, nt)

int = IntegratorFIRK(ode, getTableauSRK3(), Δt)
sol = integrate(int, nt)

pint = Integrator(pode, getTableauSymplecticEulerA(), Δt)
psol = integrate(pint, nt)

pint = Integrator(pode, getTableauSymplecticEulerB(), Δt)
psol = integrate(pint, nt)

pint = Integrator(pode, TableauIPRK(:pglrk, 2, glrk1.q, glrk1.q), Δt)
psol = integrate(pint, nt)

iint = Integrator(iode, TableauVPRK(:pglrk, 2, glrk1.q, glrk1.q, -1), Δt)
isol = integrate(iint, nt)

vint = Integrator(iode, getTableauVPLobIIIAIIIB2(), Δt)
vsol = integrate(vint, nt)

vint = Integrator(iode, getTableauVPLobIIIAIIIB3(), Δt)
vsol = integrate(vint, nt)

vint = Integrator(iode, getTableauVPLobIIIAIIIB4(), Δt)
vsol = integrate(vint, nt)

vint = IntegratorVPRKpStandard(iode, getTableauVPGLRK(1), Δt)
isol = integrate(vint, nt)

vint = IntegratorVPRKpSymplectic(iode, getTableauVPGLRK(1), Δt)
isol = integrate(vint, nt)

vint = IntegratorVPRKpSymmetric(iode, getTableauVPGLRK(1), Δt)
isol = integrate(vint, nt)

vint = IntegratorVPRKpMidpoint(iode, getTableauVPGLRK(1), Δt)
isol = integrate(vint, nt)


# TODO Add PDAE/PARK test.

dint = Integrator(idae, getTableauSymplecticProjection(:pglrk2p, glrk2.q, glrk2.q), Δt)
dsol = integrate(dint, nt)

dint = Integrator(idae, getTableauGLRKpSymplectic(2), Δt)
dsol = integrate(dint, nt)

dint = Integrator(idae, getTableauLobIIIAIIIB2pSymplectic(), Δt)
dsol = integrate(dint, nt)

dint = Integrator(idae, getTableauLobIIIAIIIB3pSymplectic(), Δt)
dsol = integrate(dint, nt)

dint = Integrator(idae, getTableauSymmetricProjection(:pglrk2p, glrk2.q, glrk2.q), Δt)
dsol = integrate(dint, nt)

dint = Integrator(idae, getTableauGLRKpSymmetric(2), Δt)
dsol = integrate(dint, nt)

dint = Integrator(idae, getTableauLobIIIAIIIB2pSymmetric(), Δt)
dsol = integrate(dint, nt)

dint = Integrator(idae, getTableauLobIIIAIIIB3pSymmetric(), Δt)
dsol = integrate(dint, nt)
