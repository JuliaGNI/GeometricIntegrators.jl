
ntime = 10
Δt    = .1
dim   = 2

x0    = rand(dim)
λ0    = [0.]

n1    = 5
x1    = rand(dim, n1)
λ1    = zeros(1, n1)

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


dae = DAE(fx, gx, fϕ, x0, λ0)
sol = Solution(dae, Δt, ntime)
@test typeof(sol) <: SolutionDAE
@test size(sol) == (2, dim, ntime+1)
@test indices(sol, 1) == 1:2
@test indices(sol, 2) == 1:dim
@test indices(sol, 3) == 0:ntime

for i in 1:ntime
    sol[1,1,i] = i
    sol[2,1,i] = i*eps()
end

@test sol.x[1,1,1] == sol[1,1,0]
@test sol.x[1:2,1:sol.nd,1] == sol[1:2,1:sol.nd,0]


dae = DAE(fx, gx, fϕ, x1, λ1)
sol = Solution(dae, Δt, ntime)
@test typeof(sol) <: SolutionDAE
@test size(sol) == (2, dim, ntime+1, n1)
@test indices(sol, 1) == 1:2
@test indices(sol, 2) == 1:dim
@test indices(sol, 3) == 0:ntime
@test indices(sol, 4) == 1:n1

for j in 1:n1
    for i in 1:ntime
        sol[1,1,i,j] = j*i
        sol[2,1,i,j] = j*i*eps()
    end
end

@test sol.x[1,1,1,1] == sol[1,1,0,1]
@test sol.x[1:2,1:sol.nd,1,1] == sol[1:2,1:sol.nd,0,1]


pdae = PDAE(fq, fp, gq, gp, gϕ, x0, x0, λ0)
psol = Solution(pdae, Δt, ntime)
@test typeof(psol) <: SolutionPDAE
@test size(psol) == (3, dim, ntime+1)
@test indices(psol, 1) == 1:3
@test indices(psol, 2) == 1:dim
@test indices(psol, 3) == 0:ntime

for i in 1:ntime
    psol[1,1,i] = i
    psol[2,1,i] = i*i
    psol[3,1,i] = i*eps()
end

@test psol.x[1,1,1] == psol[1,1,0]
@test psol.x[2,1,1] == psol[2,1,0]
@test psol.x[3,1,1] == psol[3,1,0]
@test psol.x[1:3,1:psol.nd,1] == psol[1:3,1:psol.nd,0]


pdae = PDAE(fq, fp, gq, gp, gϕ, x1, x1.^2, λ1)
psol = Solution(pdae, Δt, ntime)
@test typeof(psol) <: SolutionPDAE
@test size(psol) == (3, dim, ntime+1, n1)
@test indices(psol, 1) == 1:3
@test indices(psol, 2) == 1:dim
@test indices(psol, 3) == 0:ntime
@test indices(psol, 4) == 1:n1

for j in 1:n1
    for i in 1:ntime
        psol[1,1,i,j] = j*i
        psol[2,1,i,j] = j*i*i
        psol[3,1,i,j] = j*i*eps()
    end
end

@test psol.x[1,1,1,1] == psol[1,1,0,1]
@test psol.x[2,1,1,1] == psol[2,1,0,1]
@test psol.x[3,1,1,1] == psol[3,1,0,1]
@test psol.x[1:3,1:psol.nd,1,1] == psol[1:3,1:psol.nd,0,1]
