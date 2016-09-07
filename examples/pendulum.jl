
using GeomDAE
using PyPlot

const Δt    = 0.1
const ntime = 10000
const neps  = 1E-14
const nmax  = 20

x0 = [acos(0.4), 0.0]

function f(x, fx)
    fx[1] = x[2]
    fx[2] = sin(x[1])
end

ode = ODE(2, f, x0)

function run_pendulum(tableau, filename)
    int = Integrator(ode, tableau, Δt)
    sol = Solution(ode, ntime)

    solve!(int, sol)

    fig = figure(figsize=(6,6))
    plot(sol.x[1,:], sol.x[2,:])
    xlim(0, 6)
    ylim(-2,+2)
    savefig(filename)
end

run_pendulum(getTableauExplicitEuler(), "pendulum_explicit_euler.pdf")
run_pendulum(getTableauExplicitMidpoint(), "pendulum_explicit_midpoint.pdf")
run_pendulum(getTableauHeun(), "pendulum_heun.pdf")
run_pendulum(getTableauKutta(), "pendulum_kutta.pdf")
run_pendulum(getTableauERK4(), "pendulum_explicit_rk4.pdf")



function qf(x, fx)
    fx[:] = x
end

function pf(x, fx)
    fx[:] = sin(x)
end

function run_pendulum_partitioned(tableau, filename)

    q0 = [acos(0.4)]
    p0 = [0.0]

    ode = PODE(1, qf, pf, q0, p0)
    int = Integrator(ode, tableau, Δt)
    sol = Solution(ode, ntime)

    solve!(int, sol)

    fig = figure(figsize=(6,6))
    plot(sol.x[1,1,:], sol.x[1,2,:])
#    plot(sol.q[1,:], sol.p[1,:])
    xlim(0, 6)
    ylim(-2,+2)
    savefig(filename)
end

run_pendulum_partitioned(getTableauSymplecticEulerA(), "pendulum_symplectic_euler_a.pdf")
run_pendulum_partitioned(getTableauSymplecticEulerB(), "pendulum_symplectic_euler_b.pdf")
