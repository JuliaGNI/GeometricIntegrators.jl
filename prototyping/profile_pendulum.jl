# run with julia --track-allocation=user

using Profile
using GeometricIntegrators

set_config(:nls_solver, NewtonSolver)

const Δt    = 0.1
const ntime = 10000

function f(t, x, v)
    v[1] = x[2]
    v[2] = sin(x[1])
    nothing
end

function profile_pendulum(tableau, prof=false)
    x0 = [acos(0.4), 0.0]
    ode = ODE(f, x0)
    int = Integrator(ode, tableau, Δt)
    sol = Solution(ode, Δt, ntime)

    integrate!(int, sol)
    set_initial_conditions!(sol, ode)

    @time integrate!(int, sol)

    if prof
        set_initial_conditions!(sol, ode)

        Profile.clear()
        Profile.clear_malloc_data()

        @profile integrate!(int, sol)
    end
end


profile_pendulum(getTableauExplicitEuler())
profile_pendulum(getTableauExplicitMidpoint())
profile_pendulum(getTableauKutta())
profile_pendulum(getTableauERK4())

profile_pendulum(getTableauGLRK(1))
profile_pendulum(getTableauGLRK(2))
profile_pendulum(getTableauGLRK(3))
profile_pendulum(getTableauGLRK(4))

profile_pendulum(getTableauGLRK(1), true)
