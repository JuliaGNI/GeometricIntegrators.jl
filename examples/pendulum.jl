
using GeometricIntegrators
using PyPlot

const Δt    = 0.1
const ntime = 100000
const plot  = false

ode_tableaus = (
    (getTableauExplicitEuler(), "pendulum_explicit_euler.png"),
    (getTableauExplicitMidpoint(), "pendulum_explicit_midpoint.png"),
    (getTableauHeun(), "pendulum_heun.png"),
    (getTableauKutta(), "pendulum_kutta.png"),
    (getTableauERK4(), "pendulum_explicit_rk4.png"),
    (getTableauERK438(), "pendulum_explicit_rk438.png"),
    (getTableauImplicitEuler(), "pendulum_implicit_euler.png"),
    (getTableauGLRK1(), "pendulum_implicit_glrk1.png"),
    (getTableauGLRK2(), "pendulum_implicit_glrk2.png"),
    (getTableauGLRK3(), "pendulum_implicit_glrk3.png")
)

pode_tableaus = (
    (getTableauSymplecticEulerA(), "pendulum_symplectic_euler_a.png"),
    (getTableauSymplecticEulerB(), "pendulum_symplectic_euler_b.png"),
    (TableauEPRK(:PERK4, 4, getTableauERK4().q, getTableauERK4().q), "pendulum_explicit_prk4.png"),
    (TableauIPRK(:PGLRK2, 4, getTableauGLRK2().q, getTableauGLRK2().q), "pendulum_implicit_prk4.png")
)


function run_pendulum(tableau, filename=nothing)
    ode = pendulum_ode()
    int = Integrator(ode, tableau, Δt)
    sol = Solution(ode, Δt, ntime)

    integrate!(int, sol)
    set_initial_conditions!(sol, ode)

    print("Running ", tableau.name, "...")
    @time integrate!(int, sol)

    if filename ≠ nothing
        fig = figure(figsize=(6,6))
        plot(sol.x[1,:,1], sol.x[2,:,1])
        xlim(0, 6)
        ylim(-2,+2)
        savefig(filename)
    end
end

for tab in ode_tableaus
    if plot
        run_pendulum(tab[1], tab[2])
    else
        run_pendulum(tab[1])
    end
end


function run_pendulum_partitioned(tableau, filename=nothing)
    ode = pendulum_pode()
    int = Integrator(ode, tableau, Δt)
    sol = Solution(ode, Δt, ntime)

    integrate!(int, sol)
    set_initial_conditions!(sol, ode)

    print("Running ", tableau.name, "...")
    @time integrate!(int, sol)

    if filename ≠ nothing
        fig = figure(figsize=(6,6))
        plot(sol.x[1,1,:,1], sol.x[2,1,:,1])
        xlim(0, 6)
        ylim(-2,+2)
        savefig(filename)
    end
end

for ptab in pode_tableaus
    if plot
        run_pendulum_partitioned(ptab[1], ptab[2])
    else
        run_pendulum_partitioned(ptab[1])
    end
end

# TODO Add VPRK and VPARK examples.
