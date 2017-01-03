
Δt = 0.1
nt = 10

ode  = pendulum_ode()
pode = pendulum_pode()
iode = pendulum_iode()
idae = pendulum_idae()

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

itab = getTableauGLRK1()
iint = Integrator(iode, TableauVPRK(:pglrk, 2, itab.q, itab.q), Δt)
isol = integrate(iint, nt)

vint = Integrator(iode, getTableauLobIIIAB2(), Δt)
vsol = integrate(vint, nt)

dtab = getTableauGLRK2()
dint = Integrator(idae, getTableauSymplecticProjection(dtab.q, dtab.q), Δt)
dsol = integrate(dint, nt)
