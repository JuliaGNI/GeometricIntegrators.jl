
Δt = 0.1
ntime = 10

dim = 1
ode = ODE(dim, x -> x, [1.])
sol = Solution(ode, Δt, ntime)
@test typeof(sol) <: SolutionODE
@test size(sol) == (dim, ntime+1)
@test indices(sol, 1) == 1:dim
@test indices(sol, 2) == 0:ntime

for i in 1:ntime
    sol[1,i] = i
end

@test sol.x[1,1] == sol[1,0]
@test sol.x[1:sol.d,1] == sol[1:sol.d,0]

pode = PODE(1, (x, y) -> x, (x, y) -> 2y, [1.], [1.])
@test typeof(Solution(pode, Δt, ntime)) <: SolutionPODE

dae = DAE(2, 1, x -> x, (x, λ) -> [λ -λ], x -> x[2]-x[1], [1.; 1.], [0.])
@test typeof(Solution(dae, Δt, ntime)) <: SolutionDAE

pdae = PDAE(1, 1, (x, y) -> x, (x, y) -> 2y, (x, y, λ) -> λ, (x, y, λ) -> -λ, (x, y) -> x[1]-y[1], [1.], [1.], [0.])
@test typeof(Solution(pdae, Δt, ntime)) <: SolutionPDAE
