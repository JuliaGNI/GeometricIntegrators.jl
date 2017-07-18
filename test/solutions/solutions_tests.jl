
ntime = 10
Δt    = .1

ode   = oscillator_ode()
dim   = ode.d

q0    = rand(dim)
p0    = q0.^2
λ0    = [0.]

n1    = 5
q1    = rand(dim, n1)
p1    = q1.^2
λ1    = zeros(1, n1)


sol = Solution(ode, Δt, ntime)
@test typeof(sol) <: SolutionODE

sol0 = Solution(similar(ode, q0), Δt, ntime)
@test typeof(sol0) <: SolutionODE

sol1 = Solution(similar(ode, q1), Δt, ntime)
@test typeof(sol1) <: SolutionODE


pode = PODE(fq, fp, q0, p0)
psol = Solution(pode, Δt, ntime)

pode0 = PODE(fq, fp, q0, p0)
psol0 = Solution(pode0, Δt, ntime)
@test typeof(psol0) <: SolutionPODE

pode1 = PODE(fq, fp, q1, p1)
psol1 = Solution(pode1, Δt, ntime)
@test typeof(psol1) <: SolutionPODE


dae = DAE(fx, gx, fϕ, q0, λ0)
sol = Solution(dae, Δt, ntime)
@test typeof(sol) <: SolutionDAE

dae = DAE(fx, gx, fϕ, q1, λ1)
sol = Solution(dae, Δt, ntime)
@test typeof(sol) <: SolutionDAE


pdae = PDAE(fq, fp, gq, gp, gϕ, q0, p0, λ0)
psol = Solution(pdae, Δt, ntime)
@test typeof(psol) <: SSolutionPDAE

pdae = PDAE(fq, fp, gq, gp, gϕ, q1, p1, λ1)
psol = Solution(pdae, Δt, ntime)
@test typeof(psol) <: SSolutionPDAE
