
Δt = 0.1

@test typeof(Integrator(ODE(fx, [1.]), getTableauExplicitMidpoint(), Δt)) <: IntegratorERK
@test typeof(Integrator(ODE(fx, [1.]), getTableauCrouzeix(), Δt)) <: IntegratorDIRK
@test typeof(Integrator(ODE(fx, [1.]), getTableauImplicitMidpoint(), Δt)) <: IntegratorFIRK

ode = ODE(fx, [1.])
int = Integrator(ode, getTableauERK4(), Δt)
sol = integrate(int, 10)


pode = PODE(fq, fp, [1.], [1.])
pint = Integrator(pode, getTableauSymplecticEulerA(), Δt)
psol = integrate(int, 10)
