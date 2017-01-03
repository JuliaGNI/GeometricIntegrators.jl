
using GeometricIntegrators
using PyPlot

const Δt    = 0.1
const ntime = 100000
const splot = true


erk4  = getTableauERK4()
glrk2 = getTableauGLRK2()

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
    (TableauEPRK(:PERK4, erk4.o, erk4.q, erk4.q), "pendulum_explicit_prk4.png"),
    (TableauIPRK(:PGLRK2, glrk2.o, glrk2.q, glrk2.q), "pendulum_implicit_prk4.png")
)

iode_tableaus = (
    (TableauVPRK(:pglrk2, glrk2.o, glrk2.q, glrk2.q), "pendulum_vprk_PGLRK2.png"),
    (getTableauLobIIIAB2(), "pendulum_vprk_LobIIIAB2.png")
)

idae_tableaus = (
    (getTableauSymplecticProjection(:pglrk2p, glrk2.q, glrk2.q), "pendulum_vpark_PGLRK2.png"),
    (getTableauLobIIIAB2p(), "pendulum_vpark_LobIIIAB2.png")
)


function run_pendulum(tableau, nt; ode = pendulum_ode(), filename=nothing)
    int = Integrator(ode, tableau, Δt)
    sol = Solution(ode, Δt, nt)

    integrate!(int, sol)
    set_initial_conditions!(sol, ode)

    print("Running ", tableau.name, "...")
    @time integrate!(int, sol)

    if filename ≠ nothing
        fig = figure(figsize=(6,6))
        plot(sol.q.d[1,:], sol.q.d[2,:])
        xlim(0, 6)
        ylim(-2,+2)
        savefig(filename)
    end
end

for tab in ode_tableaus
    if splot
        run_pendulum(tab[1], ntime, filename=tab[2])
    else
        run_pendulum(tab[1], ntime)
    end
end


function run_pendulum_partitioned(tableau, nt; ode=pendulum_pode(), filename=nothing)
    int = Integrator(ode, tableau, Δt)
    sol = Solution(ode, Δt, nt)

    integrate!(int, sol)
    set_initial_conditions!(sol, ode)

    print("Running ", tableau.name, "...")
    @time integrate!(int, sol)

    if filename ≠ nothing
        fig = figure(figsize=(6,6))
        plot(sol.q.d[1,:], sol.p.d[1,:])
        xlim(0, 6)
        ylim(-2,+2)
        savefig(filename)
    end
end

for ptab in pode_tableaus
    if splot
        run_pendulum_partitioned(ptab[1], ntime, filename=ptab[2])
    else
        run_pendulum_partitioned(ptab[1], ntime)
    end
end

for ptab in iode_tableaus
    if splot
        run_pendulum_partitioned(ptab[1], 1000; ode=pendulum_iode(), filename=ptab[2])
    else
        run_pendulum_partitioned(ptab[1], 1000; ode=pendulum_iode())
    end
end

for ptab in idae_tableaus
    if splot
        run_pendulum_partitioned(ptab[1], ntime; ode=pendulum_idae(), filename=ptab[2])
    else
        run_pendulum_partitioned(ptab[1], ntime; ode=pendulum_idae())
    end
end
