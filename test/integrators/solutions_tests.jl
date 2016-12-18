
ntime = 10
Δt    = .1
dim   = 1


x0    = [1.]

ode = ODE(fx, x0)
sol = Solution(ode, Δt, ntime)
@test typeof(sol) <: SolutionODE
@test size(sol) == (dim, ntime+1)
@test indices(sol, 1) == 1:dim
@test indices(sol, 2) == 0:ntime

for i in 1:ntime
    sol[1,i] = i
end

@test sol.x[1,1] == sol[1,0]
@test sol.x[1:sol.nd,1] == sol[1:sol.nd,0]


n1    = 5
x1    = rand(1, n1)

ode = ODE(fx, x1)
sol = Solution(ode, Δt, ntime)
@test typeof(sol) <: SolutionODE
@test size(sol) == (dim, ntime+1, n1)
@test indices(sol, 1) == 1:dim
@test indices(sol, 2) == 0:ntime
@test indices(sol, 3) == 1:n1

for j in 1:n1
    for i in 1:ntime
        sol[1,i,j] = j*i
    end
end

@test sol.x[1,1,1] == sol[1,0,1]
@test sol.x[1:sol.nd,1,1] == sol[1:sol.nd,0,1]


@test typeof(Integrator(ODE(fx, x0), getTableauExplicitMidpoint(), Δt)) <: IntegratorERK


pode = PODE(fq, fp, x0, x0)
psol = Solution(pode, Δt, ntime)
@test typeof(psol) <: SolutionPODE
@test size(psol) == (2, dim, ntime+1)
@test indices(psol, 1) == 1:2
@test indices(psol, 2) == 1:dim
@test indices(psol, 3) == 0:ntime

for i in 1:ntime
    psol[1,1,i] = i
    psol[2,1,i] = i*i
end

@test psol.x[1,1,1] == psol[1,1,0]
@test psol.x[2,1,1] == psol[2,1,0]
@test psol.x[1:2,1:psol.nd,1] == psol[1:2,1:psol.nd,0]


pode = PODE(fq, fp, x1, x1)
psol = Solution(pode, Δt, ntime)
@test typeof(psol) <: SolutionPODE
@test size(psol) == (2, dim, ntime+1, n1)
@test indices(psol, 1) == 1:2
@test indices(psol, 2) == 1:dim
@test indices(psol, 3) == 0:ntime
@test indices(psol, 4) == 1:n1

for j in 1:n1
    for i in 1:ntime
        psol[1,1,i,j] = j*i
        psol[2,1,i,j] = j*i*i
    end
end

@test psol.x[1,1,1,1] == psol[1,1,0,1]
@test psol.x[2,1,1,1] == psol[2,1,0,1]
@test psol.x[1:2,1:psol.nd,1,1] == psol[1:2,1:psol.nd,0,1]


dae = DAE(x -> x, (x, λ) -> [λ -λ], x -> x[2]-x[1], [1.; 1.], [0.])
@test typeof(Solution(dae, Δt, ntime)) <: SolutionDAE

pdae = PDAE(y -> y, x -> 2x, (x, y, λ) -> λ, (x, y, λ) -> -λ, (x, y) -> x[1]-y[1], [1.], [1.], [0.])
@test typeof(Solution(pdae, Δt, ntime)) <: SolutionPDAE
