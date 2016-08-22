
@test typeof(Integrator(ODE(1, x -> x, [1.]), getTableauExplicitMidpoint())) <: IntegratorERK
@test typeof(Integrator(ODE(1, x -> x, [1.]), getTableauCrouzeix())) <: IntegratorIRK
@test typeof(Integrator(ODE(1, x -> x, [1.]), getTableauImplicitMidpoint())) <: IntegratorNLIRK

ode = ODE(1, x -> x, [1.])
int = Integrator(ode, getTableauERK4())
sol = solve(int, 10)
