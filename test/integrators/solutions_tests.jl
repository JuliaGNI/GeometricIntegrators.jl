
Δt = 0.1
ntime = 10

dim = 1
ode = ODE(dim, x -> x, [1.])
sol = Solution(ode, Δt, ntime)
@test typeof(sol) <: SolutionODE

pode = PODE(1, (x, y) -> x, (x, y) -> 2y, [1.], [1.])
@test typeof(Solution(pode, Δt, ntime)) <: SolutionPODE

dae = DAE(2, 1, x -> x, (x, λ) -> [λ -λ], x -> x[2]-x[1], [1.; 1.], [0.])
@test typeof(Solution(dae, Δt, ntime)) <: SolutionDAE

pdae = PDAE(1, 1, (x, y) -> x, (x, y) -> 2y, (x, y, λ) -> λ, (x, y, λ) -> -λ, (x, y) -> x[1]-y[1], [1.], [1.], [0.])
@test typeof(Solution(pdae, Δt, ntime)) <: SolutionPDAE
