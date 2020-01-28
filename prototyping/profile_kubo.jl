# for profiling run with julia --track-allocation=user

using Profile
using GeometricIntegrators
using GeometricIntegrators.TestProblems.KuboOscillatorProblem

sde1   = kubo_oscillator_sde_1()
sde2   = kubo_oscillator_sde_2()
sde3   = kubo_oscillator_sde_3()

using GeometricIntegrators.TestProblems.KuboOscillatorProblem: Δt

const ntime = 10000


function profile_kubo(tableau, sde, prof=false)
    int = Integrator(sde, tableau, Δt)
    sol = Solution(sde, Δt, ntime)

    integrate!(int, sol)
    set_initial_conditions!(sol, sde)

    @time integrate!(int, sol)

    if prof
        set_initial_conditions!(sol, sde)

        Profile.clear()
        Profile.clear_malloc_data()

        @profile integrate!(int, sol)
    end
end


profile_kubo(getTableauBurrageE1(), kubo_oscillator_sde_1())
profile_kubo(getTableauBurrageE1(), kubo_oscillator_sde_2())
profile_kubo(getTableauBurrageE1(), kubo_oscillator_sde_3())

profile_kubo(getTableauStochasticGLRK(1), kubo_oscillator_sde_1())
profile_kubo(getTableauStochasticGLRK(1), kubo_oscillator_sde_2())
profile_kubo(getTableauStochasticGLRK(1), kubo_oscillator_sde_3())

profile_kubo(getTableauRosslerRS1(), kubo_oscillator_sde_1())
profile_kubo(getTableauRosslerRS1(), kubo_oscillator_sde_2())
profile_kubo(getTableauRosslerRS1(), kubo_oscillator_sde_3())

profile_kubo(getTableauSRKw1(), kubo_oscillator_sde_1())
profile_kubo(getTableauSRKw1(), kubo_oscillator_sde_2())
profile_kubo(getTableauSRKw1(), kubo_oscillator_sde_3())


# profile_kubo(getTableauBurrageE1(), kubo_oscillator_sde_1(), true)
# profile_kubo(getTableauBurrageE1(), kubo_oscillator_sde_2(), true)
# profile_kubo(getTableauBurrageE1(), kubo_oscillator_sde_3(), true)
#
# profile_kubo(getTableauStochasticGLRK(1), kubo_oscillator_sde_1(), true)
# profile_kubo(getTableauStochasticGLRK(1), kubo_oscillator_sde_2(), true)
# profile_kubo(getTableauStochasticGLRK(1), kubo_oscillator_sde_3(), true)
#
# profile_kubo(getTableauRosslerRS1(), kubo_oscillator_sde_1(), true)
# profile_kubo(getTableauRosslerRS1(), kubo_oscillator_sde_2(), true)
# profile_kubo(getTableauRosslerRS1(), kubo_oscillator_sde_3(), true)
#
# profile_kubo(getTableauSRKw1(), kubo_oscillator_sde_1(), true)
# profile_kubo(getTableauSRKw1(), kubo_oscillator_sde_2(), true)
# profile_kubo(getTableauSRKw1(), kubo_oscillator_sde_3(), true)
