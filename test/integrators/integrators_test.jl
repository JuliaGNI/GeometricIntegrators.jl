
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

pint = Integrator(pode, getTableauSymplecticEulerA(), Δt)
psol = integrate(pint, nt)

pint = Integrator(pode, TableauIPRK(:pglrk, 2, glrk1.q, glrk1.q), Δt)
psol = integrate(pint, nt)

iint = Integrator(iode, TableauVPRK(:pglrk, 2, glrk1.q, glrk1.q), Δt)
isol = integrate(iint, nt)

vint = Integrator(iode, getTableauLobIIIAB2(), Δt)
vsol = integrate(vint, nt)

dint = Integrator(idae, getTableauSymplecticProjection(:pglrk2p, glrk2.q, glrk2.q), Δt)
dsol = integrate(dint, nt)

dint = Integrator(idae, getTableauLobIIIAB2p(), Δt)
dsol = integrate(dint, nt)
