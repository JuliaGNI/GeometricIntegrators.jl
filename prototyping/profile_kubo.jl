# for profiling run with julia --track-allocation=user

using Profile
using GeometricIntegrators
using GeometricIntegrators.TestProblems.KuboOscillatorProblem
using GeometricIntegrators.TestProblems.KuboOscillatorProblem: Δt

const ntime = 10000

set_config(:nls_stol_break, 1E3)


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


# SERK
profile_kubo(getTableauBurrageE1(), kubo_oscillator_sde_1())
profile_kubo(getTableauBurrageE1(), kubo_oscillator_sde_2())
profile_kubo(getTableauBurrageE1(), kubo_oscillator_sde_3())
#
# SIRK
profile_kubo(getTableauStochasticGLRK(1), kubo_oscillator_sde_1())
profile_kubo(getTableauStochasticGLRK(1), kubo_oscillator_sde_2())
profile_kubo(getTableauStochasticGLRK(1), kubo_oscillator_sde_3())
#
# WERK
profile_kubo(getTableauRosslerRS1(), kubo_oscillator_sde_1())
profile_kubo(getTableauRosslerRS1(), kubo_oscillator_sde_2())
profile_kubo(getTableauRosslerRS1(), kubo_oscillator_sde_3())

# WIRK
profile_kubo(getTableauSRKw1(), kubo_oscillator_sde_1())
profile_kubo(getTableauSRKw1(), kubo_oscillator_sde_2())
profile_kubo(getTableauSRKw1(), kubo_oscillator_sde_3())

# SIPRK
profile_kubo(getTableauStochasticStormerVerlet(), kubo_oscillator_psde_1())
profile_kubo(getTableauStochasticStormerVerlet(), kubo_oscillator_psde_2())
profile_kubo(getTableauStochasticStormerVerlet(), kubo_oscillator_psde_3())

# SISPRK
profile_kubo(getTableauModifiedStochasticStormerVerlet(), kubo_oscillator_spsde_1())
profile_kubo(getTableauModifiedStochasticStormerVerlet(), kubo_oscillator_spsde_2())
profile_kubo(getTableauModifiedStochasticStormerVerlet(), kubo_oscillator_spsde_3())



# SERK
# profile_kubo(getTableauBurrageE1(), kubo_oscillator_sde_1(), true)
# profile_kubo(getTableauBurrageE1(), kubo_oscillator_sde_2(), true)
# profile_kubo(getTableauBurrageE1(), kubo_oscillator_sde_3(), true)

# SIRK
# profile_kubo(getTableauStochasticGLRK(1), kubo_oscillator_sde_1(), true)
# profile_kubo(getTableauStochasticGLRK(1), kubo_oscillator_sde_2(), true)
# profile_kubo(getTableauStochasticGLRK(1), kubo_oscillator_sde_3(), true)

# WERK
# profile_kubo(getTableauRosslerRS1(), kubo_oscillator_sde_1(), true)
# profile_kubo(getTableauRosslerRS1(), kubo_oscillator_sde_2(), true)
# profile_kubo(getTableauRosslerRS1(), kubo_oscillator_sde_3(), true)

# WIRK
# profile_kubo(getTableauSRKw1(), kubo_oscillator_sde_1(), true)
# profile_kubo(getTableauSRKw1(), kubo_oscillator_sde_2(), true)
# profile_kubo(getTableauSRKw1(), kubo_oscillator_sde_3(), true)

# SIPRK
# profile_kubo(getTableauStochasticStormerVerlet(), kubo_oscillator_psde_1(), true)
# profile_kubo(getTableauStochasticStormerVerlet(), kubo_oscillator_psde_2(), true)
# profile_kubo(getTableauStochasticStormerVerlet(), kubo_oscillator_psde_3(), true)

# SISPRK
# profile_kubo(getTableauModifiedStochasticStormerVerlet(), kubo_oscillator_spsde_1(), true)
# profile_kubo(getTableauModifiedStochasticStormerVerlet(), kubo_oscillator_spsde_2(), true)
# profile_kubo(getTableauModifiedStochasticStormerVerlet(), kubo_oscillator_spsde_3(), true)
