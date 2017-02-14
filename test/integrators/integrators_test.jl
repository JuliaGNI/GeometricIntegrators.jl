
Δt = 0.1
nt = 10

ode  = pendulum_ode()
pode = pendulum_pode()
iode = pendulum_iode()
idae = pendulum_idae()

glrk1 = getTableauGLRK1()
glrk2 = getTableauGLRK2()

@test typeof(Integrator(ode, getTableauExplicitMidpoint(), Δt)) <: IntegratorERK
@test typeof(Integrator(ode, getTableauCrouzeix(), Δt)) <: IntegratorDIRK
@test typeof(Integrator(ode, getTableauImplicitMidpoint(), Δt)) <: IntegratorFIRK

int = Integrator(ode, getTableauExplicitEuler(), Δt)
sol = integrate(int, nt)

int = Integrator(ode, getTableauImplicitEuler(), Δt)
sol = integrate(int, nt)

int = IntegratorFIRK(ode, getTableauImplicitMidpoint(), Δt, nonlinear_solver=NewtonSolver)
sol = integrate(int, nt)

int = IntegratorFIRK(ode, getTableauGLRK1(), Δt, nonlinear_solver=NewtonSolver)
sol = integrate(int, nt)

int = IntegratorFIRK(ode, getTableauGLRK2(), Δt, nonlinear_solver=NewtonSolver)
sol = integrate(int, nt)

int = IntegratorFIRK(ode, getTableauGLRK3(), Δt, nonlinear_solver=NewtonSolver)
sol = integrate(int, nt)

int = IntegratorFIRK(ode, getTableauGLRK(4), Δt, nonlinear_solver=NewtonSolver)
sol = integrate(int, nt)

int = IntegratorFIRK(ode, getTableauSRK3(), Δt, nonlinear_solver=NewtonSolver)
sol = integrate(int, nt)

pint = Integrator(pode, getTableauSymplecticEulerA(), Δt)
psol = integrate(pint, nt)

pint = Integrator(pode, TableauIPRK(:pglrk, 2, glrk1.q, glrk1.q), Δt)
psol = integrate(pint, nt)

iint = Integrator(iode, TableauVPRK(:pglrk, 2, glrk1.q, glrk1.q, -1), Δt)
isol = integrate(iint, nt)

vint = Integrator(iode, getTableauLobIIIAIIIB2(), Δt)
vsol = integrate(vint, nt)

vint = Integrator(iode, getTableauLobIIIAIIIB3(), Δt)
vsol = integrate(vint, nt)

vint = IntegratorVPRKpStandard(iode, getTableauVPGLRK(1), Δt)
isol = integrate(vint, nt)

vint = IntegratorVPRKpSymplectic(iode, getTableauVPGLRK(1), Δt)
isol = integrate(vint, nt)

vint = IntegratorVPRKpSymmetric(iode, getTableauVPGLRK(1), Δt)
isol = integrate(vint, nt)

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
