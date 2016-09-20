
Δt = 0.1

function f(x, fx)
    fx[:] = x
end

function g(x, gx)
    gx[:] = x.^2
end

@test typeof(Integrator(ODE(1, f, [1.]), getTableauExplicitMidpoint(), Δt)) <: IntegratorERK
@test typeof(Integrator(ODE(1, f, [1.]), getTableauCrouzeix(), Δt)) <: IntegratorDIRK
@test typeof(Integrator(ODE(1, f, [1.]), getTableauImplicitMidpoint(), Δt)) <: IntegratorFIRK

ode = ODE(1, f, [1.])
int = Integrator(ode, getTableauERK4(), Δt)
sol = integrate(int, 10)


pode = PODE(1, f, g, [1.], [1.])
pint = Integrator(pode, getTableauSymplecticEulerA(), Δt)
psol = integrate(int, 10)
