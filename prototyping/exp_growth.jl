# load library
using GeometricIntegrators

# define vector field
f = (t, x, ẋ) -> ẋ[1] = x[1]

# create ODE for x = f(t,x) with initial condition x(0) = 1.0
ode = ODE(f, [1.0])

# create integrator for explicit Euler method with time step 0.1
int = Integrator(ode, getTableauExplicitEuler(), 0.1)

# compute 10 time steps and retrieve solution
sol = integrate(int, 10)

# plot and compare with the exact solution
using Plots
plot(xlim=[0,1], xlab="t", ylab="x(t)", legend=:bottomright)
plot!(sol.t.t, sol.q.d[1,:], label="numeric")
plot!(sol.t.t, exp.(sol.t.t), label="exact")

savefig("exponential_growth.pdf")
