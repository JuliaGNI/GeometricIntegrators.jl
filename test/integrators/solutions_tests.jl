
ntime = 10

x0  = [1.]
dim = 1

ode = ODE{eltype(x0)}(dim, 1, fx, 0, x0)
sol = Solution(ode, Δt, ntime)
@test typeof(sol) <: SolutionODE
# @test length(sol) == dim*(ntime+1)
@test size(sol) == (dim, ntime+1, 1)
@test indices(sol, 1) == 1:dim
@test indices(sol, 2) == 0:ntime
@test indices(sol, 3) == 1:1

for i in 1:ntime
    sol[1,i] = i
end

@test sol.x[1,1,1] == sol[1,0,1]
@test sol.x[1:sol.nd,1,1] == sol[1:sol.nd,0,1]

@test typeof(Integrator(ODE(fx, [1.]), getTableauExplicitMidpoint(), Δt)) <: IntegratorERK


pode = PODE(fq, fp, [1.], [1.])
psol = Solution(pode, Δt, ntime)
@test typeof(psol) <: SolutionPODE
@test size(psol) == (2, dim, ntime+1, 1)
@test indices(psol, 1) == 1:2
@test indices(psol, 2) == 1:dim
@test indices(psol, 3) == 0:ntime
@test indices(psol, 4) == 1:1

for i in 1:ntime
    psol[1,1,i,1] = i
    psol[2,1,i,1] = i*i
end

@test psol.x[1,1,1,1] == psol[1,1,0,1]
@test psol.x[2,1,1,1] == psol[2,1,0,1]
@test psol.x[1:2,1:psol.nd,1,1] == psol[1:2,1:psol.nd,0,1]


dae = DAE(x -> x, (x, λ) -> [λ -λ], x -> x[2]-x[1], [1.; 1.], [0.])
@test typeof(Solution(dae, Δt, ntime)) <: SolutionDAE

pdae = PDAE(y -> y, x -> 2x, (x, y, λ) -> λ, (x, y, λ) -> -λ, (x, y) -> x[1]-y[1], [1.], [1.], [0.])
@test typeof(Solution(pdae, Δt, ntime)) <: SolutionPDAE
