
ntime = 10
Δt    = .1
dim   = 2

q0    = rand(dim)
p0    = q0.^2
λ0    = [0.]

n1    = 5
q1    = rand(dim, n1)
p1    = q1.^2
λ1    = zeros(1, n1)

ode = ODE(fx, q0)
sol = Solution(ode, Δt, ntime)
@test typeof(sol) <: SolutionODE

ode = ODE(fx, q1)
sol = Solution(ode, Δt, ntime)
@test typeof(sol) <: SolutionODE

@test typeof(Integrator(ODE(fx, q0), getTableauExplicitMidpoint(), Δt)) <: IntegratorERK


pode = PODE(fq, fp, q0, p0)
psol = Solution(pode, Δt, ntime)
@test typeof(psol) <: SolutionPODE

pode = PODE(fq, fp, q1, p1)
psol = Solution(pode, Δt, ntime)
@test typeof(psol) <: SolutionPODE


dae = DAE(fx, gx, fϕ, q0, λ0)
sol = Solution(dae, Δt, ntime)
@test typeof(sol) <: SolutionDAE

dae = DAE(fx, gx, fϕ, q1, λ1)
sol = Solution(dae, Δt, ntime)
@test typeof(sol) <: SolutionDAE


pdae = PDAE(fq, fp, gq, gp, gϕ, q0, p0, λ0)
psol = Solution(pdae, Δt, ntime)
@test typeof(psol) <: SolutionPDAE

pdae = PDAE(fq, fp, gq, gp, gϕ, q1, p1, λ1)
psol = Solution(pdae, Δt, ntime)
@test typeof(psol) <: SolutionPDAE
