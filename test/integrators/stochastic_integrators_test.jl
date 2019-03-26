
set_config(:nls_solver, NewtonSolver)


function rel_energy_err_sde(sol)

    NQ = ndims(sol.q)

    if NQ==2

        en_ref  = 0.5 * ( sol.q.d[1,1]^2 + sol.q.d[2,1]^2 )
        en_last = 0.5 * ( sol.q.d[1,end]^2 + sol.q.d[2,end]^2 )

        return abs((en_last-en_ref)/en_ref)

    elseif NQ==3

        en_ref  = 0.5 * ( sol.q.d[1,1,:].^2 + sol.q.d[2,1,:].^2 )
        en_last = 0.5 * ( sol.q.d[1,end,:].^2 + sol.q.d[2,end,:].^2 )

        return maximum(abs.( (en_last .- en_ref) ./ en_ref ))

    elseif NQ==4

        en_ref  = 0.5 * ( sol.q.d[1,1,:,:].^2 + sol.q.d[2,1,:,:].^2 )
        en_last = 0.5 * ( sol.q.d[1,end,:,:].^2 + sol.q.d[2,end,:,:].^2 )

        return maximum(abs.( (en_last .- en_ref) ./ en_ref ))

    end
end


function rel_energy_err_psde(sol)

    NQ = ndims(sol.q)

    if NQ==2

        en_ref  = 0.5 * ( sol.q.d[1,1]^2 + sol.p.d[1,1]^2 )
        en_last = 0.5 * ( sol.q.d[1,end]^2 + sol.p.d[1,end]^2 )

        return abs((en_last-en_ref)/en_ref)

    elseif NQ==3

        en_ref  = 0.5 * ( sol.q.d[1,1,:].^2 + sol.p.d[1,1,:].^2 )
        en_last = 0.5 * ( sol.q.d[1,end,:].^2 + sol.p.d[1,end,:].^2 )

        return maximum(abs.( (en_last .- en_ref) ./ en_ref ))

    elseif NQ==4

        en_ref  = 0.5 * ( sol.q.d[1,1,:,:].^2 + sol.p.d[1,1,:,:].^2 )
        en_last = 0.5 * ( sol.q.d[1,end,:,:].^2 + sol.p.d[1,end,:,:].^2 )

        return maximum(abs.( (en_last .- en_ref) ./ en_ref ))

    end
end


Δt  = 0.01
nt  = 10

sde1   = kubo_oscillator_sde_1()
sde2   = kubo_oscillator_sde_2()
sde3   = kubo_oscillator_sde_3()
psde1  = kubo_oscillator_psde_1()
psde2  = kubo_oscillator_psde_2()
psde3  = kubo_oscillator_psde_3()
spsde1  = kubo_oscillator_spsde_1()
spsde2  = kubo_oscillator_spsde_2()
spsde3  = kubo_oscillator_spsde_3()


@test typeof(Integrator(sde1, getTableauBurrageE1(), Δt)) <: IntegratorSERK
@test typeof(Integrator(sde1, getTableauStochasticGLRK(1), Δt)) <: IntegratorSIRK
@test typeof(Integrator(sde1, getTableauRosslerRS1(), Δt)) <: IntegratorWERK
@test typeof(Integrator(sde1, getTableauSRKw1(), Δt)) <: IntegratorWIRK
@test typeof(Integrator(psde1, getTableauStochasticStormerVerlet(), Δt)) <: IntegratorSIPRK
@test typeof(Integrator(spsde1, getTableauModifiedStochasticStormerVerlet(), Δt)) <: IntegratorSISPRK


### SERK Integrators ###
int = Integrator(sde1, getTableauBurrageE1(), Δt)
sol = Solution(sde1, Δt, nt, conv="strong")
integrate!(int, sol)

@test rel_energy_err_sde(sol) < 2E-6

int = Integrator(sde2, getTableauBurrageE1(), Δt)
sol = Solution(sde2, Δt, nt, conv="strong")
integrate!(int, sol)

@test rel_energy_err_sde(sol) < 2E-6

int = Integrator(sde3, getTableauBurrageE1(), Δt)
sol = Solution(sde3, Δt, nt, conv="strong")
integrate!(int, sol)

@test rel_energy_err_sde(sol) < 2E-6


### SIRK Integrators ###
int = Integrator(sde1, getTableauStochasticGLRK(1), Δt)
sol = Solution(sde1, Δt, nt, conv="strong")
integrate!(int, sol)

@test rel_energy_err_sde(sol) < 1E-14

int = Integrator(sde2, getTableauStochasticGLRK(1), Δt)
sol = Solution(sde2, Δt, nt, conv="strong")
integrate!(int, sol)

@test rel_energy_err_sde(sol) < 1E-14

int = Integrator(sde3, getTableauStochasticGLRK(1), Δt)
sol = Solution(sde3, Δt, nt, conv="strong")
integrate!(int, sol)

@test rel_energy_err_sde(sol) < 1E-14


### WERK Integrators ###
int = Integrator(sde1, getTableauRosslerRS1(), Δt)
sol = Solution(sde1, Δt, nt, conv="weak")
integrate!(int, sol)

@test rel_energy_err_sde(sol) < 1E-5

int = Integrator(sde2, getTableauRosslerRS1(), Δt)
sol = Solution(sde2, Δt, nt, conv="weak")
integrate!(int, sol)

@test rel_energy_err_sde(sol) < 1E-5

int = Integrator(sde3, getTableauRosslerRS1(), Δt)
sol = Solution(sde3, Δt, nt, conv="weak")
integrate!(int, sol)

@test rel_energy_err_sde(sol) < 1E-5

### WIRK Integrators ###
int = Integrator(sde1, getTableauSRKw1(), Δt)
sol = Solution(sde1, Δt, nt, conv="weak")
integrate!(int, sol)

@test rel_energy_err_sde(sol) < 1E-14

int = Integrator(sde2, getTableauSRKw1(), Δt)
sol = Solution(sde2, Δt, nt, conv="weak")
integrate!(int, sol)

@test rel_energy_err_sde(sol) < 1E-14

int = Integrator(sde3, getTableauSRKw1(), Δt)
sol = Solution(sde3, Δt, nt, conv="weak")
integrate!(int, sol)

@test rel_energy_err_sde(sol) < 1E-14


### SIPRK Integrators ###
int = Integrator(psde1, getTableauStochasticStormerVerlet(), Δt)
sol = Solution(psde1, Δt, nt, conv="strong")
integrate!(int, sol)

@test rel_energy_err_psde(sol) < 1E-5

int = Integrator(psde2, getTableauStochasticStormerVerlet(), Δt)
sol = Solution(psde2, Δt, nt, conv="strong")
integrate!(int, sol)

@test rel_energy_err_psde(sol) < 1E-5

int = Integrator(psde3, getTableauStochasticStormerVerlet(), Δt)
sol = Solution(psde3, Δt, nt, conv="strong")
integrate!(int, sol)

@test rel_energy_err_psde(sol) < 1E-5


### SISPRK Integrators ###
int = Integrator(spsde1, getTableauModifiedStochasticStormerVerlet(), Δt)
sol = Solution(spsde1, Δt, nt, conv="strong")
integrate!(int, sol)

@test rel_energy_err_psde(sol) < 0.02

int = Integrator(spsde2, getTableauModifiedStochasticStormerVerlet(), Δt)
sol = Solution(spsde2, Δt, nt, conv="strong")
integrate!(int, sol)

@test rel_energy_err_psde(sol) < 0.02

int = Integrator(spsde3, getTableauModifiedStochasticStormerVerlet(), Δt)
sol = Solution(spsde3, Δt, nt, conv="strong")
integrate!(int, sol)

@test rel_energy_err_psde(sol) < 0.02
