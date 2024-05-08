# using GeometricSolutions: test_interface
using GeometricProblems.HarmonicOscillator
using GeometricIntegrators
using Test


nt = 10
Δt = 0.1
ni = 5

t0 = 0.0
x0 = rand(2)
y0 = rand(2)
q0 = rand(1)
p0 = rand(1)
z0 = rand(3)
λ0 = rand(1)
μ0 = rand(1)

t1 = 1.0
x1 = [rand(2) for j = 1:ni]
y1 = [rand(2) for j = 1:ni]
q1 = [rand(1) for j = 1:ni]
p1 = [rand(1) for j = 1:ni]
z1 = [rand(3) for j = 1:ni]
λ1 = [rand(1) for j = 1:ni]
μ1 = [rand(1) for j = 1:ni]

t2 = t1 + (t1 - t0)
x2 = [rand(2) for j = 1:ni]
y2 = [rand(2) for j = 1:ni]
q2 = [rand(1) for j = 1:ni]
p2 = [rand(1) for j = 1:ni]
z2 = [rand(3) for j = 1:ni]
λ2 = [rand(1) for j = 1:ni]
μ2 = [rand(1) for j = 1:ni]

tx = zero(x0)
ty = zero(y0)
tq = zero(q0)
tp = zero(p0)
tz = zero(z0)
tλ = zero(λ0)
tμ = zero(μ0)

xs = [rand(2) for i = 1:nt]
ys = [rand(2) for i = 1:nt]
qs = [rand(1) for i = 1:nt]
ps = [rand(1) for i = 1:nt]
zs = [rand(3) for i = 1:nt]
λs = [rand(1) for i = 1:nt]
μs = [rand(1) for i = 1:nt]

Xs = [rand(2) for i = 1:nt, j = 1:ni]
Ys = [rand(2) for i = 1:nt, j = 1:ni]
Qs = [rand(1) for i = 1:nt, j = 1:ni]
Ps = [rand(1) for i = 1:nt, j = 1:ni]
Zs = [rand(3) for i = 1:nt, j = 1:ni]
Λs = [rand(1) for i = 1:nt, j = 1:ni]
Ms = [rand(1) for i = 1:nt, j = 1:ni]

ode  = odeproblem()
dae  = daeproblem()
pode = podeproblem()
pdae = pdaeproblem()


@testset "$(rpad("ODE Solution",80))" begin
    solstep = SolutionStep(ode)

    # test constructors and general functionality
    sol = Solution(ode)
    @test typeof(sol) <: SolutionODE

    # test_interface(sol) # TODO reactivate

    sol0 = Solution(similar(ode, ics = (q = StateVariable(x0),)))
    @test typeof(sol0) <: SolutionODE
    @test sol != sol0

    # test initial conditions
    copy!(solstep, sol0[0])

    @test solstep.t == t0
    @test solstep.q == x0

    # test set/get solution
    sol1 = Solution(similar(ode, ics = (q = StateVariable(x0),)))
    sol2 = Solution(similar(ode, ics = (q = StateVariable(x0),)))
    for i in eachindex(xs)
        solstep.q .= xs[i]
        sol1[i] = (q = copy(xs[i]),)
        sol2[i] = solstep
    end
    @test sol1.q[1:nt] == xs
    @test sol2.q[1:nt] == xs

    # test step and nstore parameters
    sol = Solution(similar(ode, tspan = 2 .* tspan(ode)); step = 2)
    @test ntime(sol) == 20
    @test nstore(sol) == 10

    sol = Solution(similar(ode, tspan = 2 .* tspan(ode)); step = 10)
    @test ntime(sol) == 20
    @test nstore(sol) == 2

end


@testset "$(rpad("PODE Solution",80))" begin
    solstep = SolutionStep(pode)

    # test constructors and general functionality
    sol = Solution(pode)
    @test typeof(sol) <: SolutionPODE

    # test_interface(sol)

    sol0 = Solution(similar(pode, ics = (q = StateVariable(q0), p = StateVariable(p0))))
    @test typeof(sol0) <: SolutionPODE
    @test sol != sol0

    # test initial conditions
    copy!(solstep, sol0[0])

    @test solstep.t == t0
    @test solstep.q == q0
    @test solstep.p == p0

    # test set/get solution
    sol1 = Solution(similar(pode, ics = (q = StateVariable(q0), p = StateVariable(p0))))
    sol2 = Solution(similar(pode, ics = (q = StateVariable(q0), p = StateVariable(p0))))
    for i = 1:nt
        solstep.q .= qs[i]
        solstep.p .= ps[i]
        sol1[i] = (q = copy(qs[i]), p = copy(ps[i]))
        sol2[i] = solstep
    end
    @test sol1.q[1:nt] == qs
    @test sol1.p[1:nt] == ps
    @test sol2.q[1:nt] == qs
    @test sol2.p[1:nt] == ps

    # test step and nstore parameters
    sol = Solution(similar(pode, tspan = 2 .* tspan(pode)), step = 2)
    @test ntime(sol) == 20
    @test nstore(sol) == 10

    sol = Solution(similar(pode, tspan = 2 .* tspan(pode)), step = 10)
    @test ntime(sol) == 20
    @test nstore(sol) == 2
end


@testset "$(rpad("DAE Solution",80))" begin
    solstep = SolutionStep(dae)

    # test constructors and general functionality
    sol = Solution(dae)
    @test typeof(sol) <: SolutionDAE

    # test_interface(sol)

    sol0 = Solution(similar(dae, ics = (q = StateVariable(x0), λ = AlgebraicVariable(λ0), μ = AlgebraicVariable(μ0))))
    @test typeof(sol0) <: SolutionDAE
    @test sol != sol0

    # test initial conditions
    copy!(solstep, sol0[0])

    @test solstep.t == t0
    @test solstep.q == x0
    @test solstep.λ == λ0

    # test set/get solution
    sol1 = Solution(similar(dae, ics = (q = StateVariable(x0), λ = AlgebraicVariable(λ0), μ = AlgebraicVariable(μ0))))
    sol2 = Solution(similar(dae, ics = (q = StateVariable(x0), λ = AlgebraicVariable(λ0), μ = AlgebraicVariable(μ0))))
    for i = 1:nt
        solstep.q .= xs[i]
        solstep.λ .= λs[i]
        sol1[i] = (q = copy(xs[i]), λ = copy(λs[i]))
        sol2[i] = solstep
    end
    @test sol1.q[1:nt] == xs
    @test sol1.λ[1:nt] == λs
    @test sol2.q[1:nt] == xs
    @test sol2.λ[1:nt] == λs

    # test step and nstore parameters
    sol = Solution(similar(dae, tspan = 2 .* tspan(dae)), step = 2)
    @test ntime(sol) == 20
    @test nstore(sol) == 10

    sol = Solution(similar(dae, tspan = 2 .* tspan(dae)), step = 10)
    @test ntime(sol) == 20
    @test nstore(sol) == 2
end


@testset "$(rpad("PDAE Solution",80))" begin
    solstep = SolutionStep(pdae)

    # test constructors and general functionality
    sol = Solution(pdae)
    @test typeof(sol) <: SolutionPDAE

    # test_interface(sol)

    sol0 = Solution(similar(pdae, ics = (q = StateVariable(q0), p = StateVariable(p0), λ = AlgebraicVariable(λ0), μ = AlgebraicVariable(μ0))))
    @test typeof(sol0) <: SolutionPDAE
    @test sol != sol0

    # test initial conditions
    copy!(solstep, sol0[0])

    @test solstep.t == t0
    @test solstep.q == q0
    @test solstep.p == p0
    @test solstep.λ == λ0

    # test set/get solution
    sol1 = Solution(similar(pdae, ics = (q = StateVariable(q0), p = StateVariable(p0), λ = AlgebraicVariable(λ0), μ = AlgebraicVariable(μ0))))
    sol2 = Solution(similar(pdae, ics = (q = StateVariable(q0), p = StateVariable(p0), λ = AlgebraicVariable(λ0), μ = AlgebraicVariable(μ0))))
    for i = 1:nt
        solstep.q .= qs[i]
        solstep.p .= ps[i]
        solstep.λ .= λs[i]
        sol1[i] = (q = copy(qs[i]), p = copy(ps[i]), λ = copy(λs[i]))
        sol2[i] = solstep
    end
    @test sol1.q[1:nt] == qs
    @test sol1.p[1:nt] == ps
    @test sol1.λ[1:nt] == λs
    @test sol2.q[1:nt] == qs
    @test sol2.p[1:nt] == ps

    # test step and nstore parameters
    sol = Solution(similar(pdae, tspan = 2 .* tspan(pdae)), step = 2)
    @test ntime(sol) == 20
    @test nstore(sol) == 10

    sol = Solution(similar(pdae, tspan = 2 .* tspan(pdae)), step = 10)
    @test ntime(sol) == 20
    @test nstore(sol) == 2

end
