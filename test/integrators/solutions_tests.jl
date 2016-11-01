
ntime = 10

x0  = [1.]
dim = 1

ode = ODE{eltype(x0)}(dim, fx, 0, x0)
sol = Solution(ode, ntime)
@test typeof(sol) <: SolutionODE
@test size(sol) == (dim, ntime+1)
@test indices(sol, 1) == 1:dim
@test indices(sol, 2) == 0:ntime

for i in 1:ntime
    sol[1,i] = i
end

@test sol.x[1,1] == sol[1,0]
@test sol.x[1:sol.d,1] == sol[1:sol.d,0]

@test typeof(Integrator(ODE(fx, [1.]), getTableauExplicitMidpoint(), Δt)) <: IntegratorERK


pode = PODE(fq, fp, [1.], [1.])
psol = Solution(pode, ntime)
@test typeof(psol) <: SolutionPODE
@test size(psol) == (dim, 2, ntime+1)
@test indices(psol, 1) == 1:dim
@test indices(psol, 2) == 1:2
@test indices(psol, 3) == 0:ntime

for i in 1:ntime
    psol[1,1,i] = i
    psol[1,2,i] = i*i
end

@test psol.x[1,1,1] == psol[1,1,0]
@test psol.x[1,2,1] == psol[1,2,0]
@test psol.x[1:sol.d,1:2,1] == psol[1:sol.d,1:2,0]


dae = DAE(x -> x, (x, λ) -> [λ -λ], x -> x[2]-x[1], [1.; 1.], [0.])
@test typeof(Solution(dae, ntime)) <: SolutionDAE

pdae = PDAE(y -> y, x -> 2x, (x, y, λ) -> λ, (x, y, λ) -> -λ, (x, y) -> x[1]-y[1], [1.], [1.], [0.])
@test typeof(Solution(pdae, ntime)) <: SolutionPDAE
