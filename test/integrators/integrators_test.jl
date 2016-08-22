
@test typeof(Integrator(ODE(1, x -> x, [1.]), getTableauExplicitMidpoint())) <: IntegratorERK
@test typeof(Integrator(ODE(1, x -> x, [1.]), getTableauCrouzeix())) <: IntegratorDIRK
@test typeof(Integrator(ODE(1, x -> x, [1.]), getTableauImplicitMidpoint())) <: IntegratorFIRK

ode = ODE(1, x -> x, [1.])
int = Integrator(ode, getTableauERK4())
sol = solve(int, 10)


pode = PODE(1, y -> y, x -> 2x, [1.], [1.])
pint = Integrator(pode, getTableauSymplecticEulerA())
psol = solve(int, 10)
