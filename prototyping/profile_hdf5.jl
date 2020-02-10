# run with julia --track-allocation=user

using Profile
using GeometricIntegrators

#set_config(:nls_solver, NewtonSolver)

const Δt    = 0.1
const ntime = 10000

function f(t, x, v)
    v[1] = x[2]
    v[2] = sin(x[1])
    nothing
end

function profile_pendulum(tableau, file, prof=false)
    x0 = [acos(0.4), 0.0]
    ode = ODE(f, x0)
    int = Integrator(ode, tableau, Δt)
    sol = Solution(ode, Δt, ntime; filename=file)

    integrate!(int, sol)
    write_to_hdf5(sol)

    @time set_initial_conditions!(sol, ode)
    integrate!(int, sol)
    @time write_to_hdf5(sol)

    if prof
        Profile.clear()
        Profile.clear_malloc_data()

        @profile set_initial_conditions!(sol, ode)
        integrate!(int, sol)
        @profile write_to_hdf5(sol)
    end

    close(sol)
    rm(file)
end


profile_pendulum(getTableauERK4(), "test.hdf5", true)
