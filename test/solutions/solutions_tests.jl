
ntime = 10
Δt    = .1

ode   = oscillator_ode()
dim   = ode.d

t0    = 0.
q0    = rand(dim)
p0    = q0.^2
λ0    = [0.]

t1    = 1.
n1    = 5
q1    = rand(dim, n1)
p1    = q1.^2
λ1    = zeros(1, n1)

tq    = zero(q0)
q2    = rand(dim, n1)


### Test SolutionODE ###

sol = Solution(ode, Δt, ntime)
@test typeof(sol) <: SolutionODE

sol0 = Solution(similar(ode, q0), Δt, ntime)
@test typeof(sol0) <: SolutionODE

sol1 = Solution(similar(ode, q1), Δt, ntime)
@test typeof(sol1) <: SolutionODE

@test sol != sol0
@test sol != sol1

set_initial_conditions!(sol, t0, q0)
get_initial_conditions!(sol, tq, 1)
@test sol != sol0
@test tq == q0

set_initial_conditions!(sol1, similar(ode, t1, q2))
get_initial_conditions!(sol1, tq, 1)
@test tq == q2[:,1]


### Test SolutionPODE ###

pode = PODE(fq, fp, q0, p0)
psol = Solution(pode, Δt, ntime)

pode0 = PODE(fq, fp, q0, p0)
psol0 = Solution(pode0, Δt, ntime)
@test typeof(psol0) <: SolutionPODE

pode1 = PODE(fq, fp, q1, p1)
psol1 = Solution(pode1, Δt, ntime)
@test typeof(psol1) <: SolutionPODE

@test pode == pode0
@test pode != pode1


### Test SolutionDAE ###

dae = DAE(fx, gx, fϕ, q0, λ0)
sol = Solution(dae, Δt, ntime)
@test typeof(sol) <: SolutionDAE

dae = DAE(fx, gx, fϕ, q1, λ1)
sol = Solution(dae, Δt, ntime)
@test typeof(sol) <: SolutionDAE


### Test SolutionPDAE ###

pdae = PDAE(fq, fp, gq, gp, gϕ, q0, p0, λ0)
psol = Solution(pdae, Δt, ntime)
@test typeof(psol) <: SSolutionPDAE

pdae = PDAE(fq, fp, gq, gp, gϕ, q1, p1, λ1)
psol = Solution(pdae, Δt, ntime)
@test typeof(psol) <: SSolutionPDAE


### Test SolutionSDE ###
sde  = kubo_oscillator_sde_1()
ssol = Solution(sde, Δt, ntime)
@test typeof(ssol) <: SolutionSDE

### Test SolutionPSDE ###
psde  = kubo_oscillator_psde_1()
ssol = Solution(psde, Δt, ntime)
@test typeof(ssol) <: SolutionPSDE

### Test SolutionSPSDE ###
spsde  = kubo_oscillator_spsde_1()
ssol = Solution(spsde, Δt, ntime)
@test typeof(ssol) <: SolutionPSDE
