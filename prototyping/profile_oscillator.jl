# run with julia --track-allocation=user

using Profile
using GeometricIntegrators
using GeometricIntegrators.TestProblems.HarmonicOscillatorProblem

set_config(:nls_stol_break, Inf)

const Δt = 0.1
const ntime = 10000


function profile_oscillator(tableau, equ, prof=false)
    int = Integrator(equ, tableau, Δt)
    sol = Solution(equ, Δt, ntime)

    integrate!(int, sol)
    set_initial_conditions!(sol, equ)

    @time integrate!(int, sol)

    if prof
        set_initial_conditions!(sol, equ)

        Profile.clear()
        Profile.clear_malloc_data()

        @profile integrate!(int, sol)
    end
end


profile_oscillator(getTableauExplicitEuler(), harmonic_oscillator_ode())
profile_oscillator(getTableauExplicitMidpoint(), harmonic_oscillator_ode())
profile_oscillator(getTableauKutta(), harmonic_oscillator_ode())
profile_oscillator(getTableauERK4(), harmonic_oscillator_ode())

profile_oscillator(getTableauCrouzeix(), harmonic_oscillator_ode())

profile_oscillator(getTableauGLRK(1), harmonic_oscillator_ode())
profile_oscillator(getTableauGLRK(2), harmonic_oscillator_ode())
profile_oscillator(getTableauGLRK(3), harmonic_oscillator_ode())
profile_oscillator(getTableauGLRK(4), harmonic_oscillator_ode())

profile_oscillator(getTableauSymplecticEulerA(), harmonic_oscillator_pode())
profile_oscillator(getTableauSymplecticEulerB(), harmonic_oscillator_pode())
profile_oscillator(getTableauLobattoIIIAIIIB2(), harmonic_oscillator_pode())
profile_oscillator(getTableauLobattoIIIBIIIA2(), harmonic_oscillator_pode())

profile_oscillator(TableauEPRK(:perk4, 4, getTableauERK4().q), harmonic_oscillator_pode())

profile_oscillator(TableauIPRK(:pglrk, 2, getCoefficientsGLRK(1)), harmonic_oscillator_pode())
profile_oscillator(TableauIPRK(:pglrk, 4, getCoefficientsGLRK(2)), harmonic_oscillator_pode())
profile_oscillator(TableauIPRK(:pglrk, 6, getCoefficientsGLRK(3)), harmonic_oscillator_pode())
profile_oscillator(TableauIPRK(:pglrk, 8, getCoefficientsGLRK(4)), harmonic_oscillator_pode())

# profile_oscillator(getTableauGLRK(1), true)
