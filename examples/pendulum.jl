
using GeomDAE
using PyPlot

const Δt    = 0.1
const ntime = 10000
const neps  = 1E-14
const nmax  = 20


function run_pendulum(tableau, filename)

    f  = x -> [x[2], sin(x[1])]
    x0 = [acos(0.4), 0.0]

    ode = ODE(2, f, x0)
    int = Integrator(ode, tableau)
    sol = Solution(ode, Δt, ntime)

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


function run_pendulum_partitioned(tableau, filename)

    f  = p -> p
    g  = q -> sin(q)
    q0 = [acos(0.4)]
    p0 = [0.0]

    ode = PODE(1, f, g, q0, p0)
    int = Integrator(ode, tableau)
    sol = Solution(ode, Δt, ntime)

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
