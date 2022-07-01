using GeometricSolutions: test_interface
using GeometricIntegrators
using GeometricIntegrators.Solutions: offset
using GeometricProblems.HarmonicOscillator
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
μ0 = rand(2)

t1 = 1.0
x1 = [rand(2) for j = 1:ni]
y1 = [rand(2) for j = 1:ni]
q1 = [rand(1) for j = 1:ni]
p1 = [rand(1) for j = 1:ni]
z1 = [rand(3) for j = 1:ni]
λ1 = [rand(1) for j = 1:ni]
μ1 = [rand(2) for j = 1:ni]

t2 = t1 + (t1 - t0)
x2 = [rand(2) for j = 1:ni]
y2 = [rand(2) for j = 1:ni]
q2 = [rand(1) for j = 1:ni]
p2 = [rand(1) for j = 1:ni]
z2 = [rand(3) for j = 1:ni]
λ2 = [rand(1) for j = 1:ni]
μ2 = [rand(2) for j = 1:ni]

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
μs = [rand(2) for i = 1:nt]

Xs = [rand(2) for i = 1:nt, j = 1:ni]
Ys = [rand(2) for i = 1:nt, j = 1:ni]
Qs = [rand(1) for i = 1:nt, j = 1:ni]
Ps = [rand(1) for i = 1:nt, j = 1:ni]
Zs = [rand(3) for i = 1:nt, j = 1:ni]
Λs = [rand(1) for i = 1:nt, j = 1:ni]
Ms = [rand(2) for i = 1:nt, j = 1:ni]

ode = harmonic_oscillator_ode()
dae = harmonic_oscillator_dae()
pode = harmonic_oscillator_pode()
pdae = harmonic_oscillator_pdae()


@testset "$(rpad("ODE Solution",80))" begin
    asol = AtomicSolution(ode)

    # test constructors and general functionality
    sol = Solution(ode)
    @test typeof(sol) <: SolutionODE

    # test_interface(sol) # TODO reactivate

    sol0 = Solution(similar(ode, ics=(q=x0,)))
    @test typeof(sol0) <: SolutionODE

    # sol1 = Solution(similar(ode, ics=(q=x1,)))
    # @test typeof(sol1) <: SolutionODE

    @test sol != sol0
    # @test sol != sol1

    @test nsave(sol) == 1
    @test ntime(sol) == nt
    @test timesteps(sol) == Δt .* collect(0:nt)

    # test initial conditions
    set_initial_conditions!(sol, t0, x0)
    get_initial_conditions!(sol, tx, 1)
    get_initial_conditions!(sol, asol, 1)

    @test tx == x0
    @test asol.t == t0
    @test asol.q == x0

    δt, δx = get_initial_conditions(sol, 1)
    @test δt == t0
    @test δx == x0

    # set_initial_conditions!(sol1, similar(ode, t1, x2))
    # for j = 1:ni
    #     get_initial_conditions!(sol1, tx, j)
    #     @test tx == x2[j]
    # end

    # test set/get solution
    sol1 = Solution(similar(ode, ics=(q=x0,)))
    sol2 = Solution(similar(ode, ics=(q=x0,)))
    for i = 1:nt
        tx .= xs[i]
        asol.q .= xs[i]
        set_solution!(sol1, tx, i)
        set_solution!(sol2, asol, i)
    end
    @test sol1.q[1:nt] == xs
    @test sol2.q[1:nt] == xs

    # sol1 = Solution(similar(ode, ics=(q=x1,)))
    # sol2 = Solution(similar(ode, ics=(q=x1,)))
    # for i = 1:nt
    #     for j = 1:ni
    #         tx .= Xs[i, j]
    #         asol.q .= Xs[i, j]
    #         set_solution!(sol1, tx, i, j)
    #         set_solution!(sol2, asol, i, j)
    #     end
    # end
    # @test sol1.q[1:nt, :] == Xs
    # @test sol2.q[1:nt, :] == Xs

    # test nsave and nwrite parameters
    sol = Solution(similar(ode, tspan = 2 .* tspan(ode)), nsave = 2)
    @test sol.nt == 10

    sol = Solution(similar(ode, tspan = 2 .* tspan(ode)), nsave = 2, nwrite = 10)
    @test sol.nt == 5

    # test reset
    sol = Solution(ode)
    reset!(sol)
    @test sol.t[0] == t1
    @test sol.t[end] == t2
    @test offset(sol) == nt
end


@testset "$(rpad("PODE Solution",80))" begin
    asol = AtomicSolution(pode)

    # test constructors and general functionality
    sol = Solution(pode)
    @test typeof(sol) <: SolutionPODE

    # test_interface(sol) # TODO reactivate

    sol0 = Solution(similar(pode, ics=(q=q0, p=p0)))
    @test typeof(sol0) <: SolutionPODE

    # sol1 = Solution(similar(pode, ics=(q=q1, p=p1)))
    # @test typeof(sol1) <: SolutionPODE

    @test sol != sol0
    # @test sol != sol1

    @test nsave(sol) == 1
    @test ntime(sol) == nt
    @test timesteps(sol) == Δt .* collect(0:nt)

    # test initial conditions
    set_initial_conditions!(sol, t0, q0, p0)
    get_initial_conditions!(sol, tq, tp, 1)
    get_initial_conditions!(sol, asol, 1)

    @test tq == q0
    @test tp == p0
    @test asol.t == t0
    @test asol.q == q0
    @test asol.p == p0

    δt, δq, δp = get_initial_conditions(sol, 1)
    @test δt == t0
    @test δq == q0
    @test δp == p0

    # set_initial_conditions!(sol1, similar(pode, t1, q2, p2))
    # for j = 1:ni
    #     get_initial_conditions!(sol1, tq, tp, j)
    #     @test tq == q2[j]
    #     @test tp == p2[j]
    # end

    # test set/get solution
    sol1 = Solution(similar(pode, ics=(q=q0, p=p0)))
    sol2 = Solution(similar(pode, ics=(q=q0, p=p0)))
    for i = 1:nt
        tq .= qs[i]
        tp .= ps[i]
        asol.q .= qs[i]
        asol.p .= ps[i]
        set_solution!(sol1, tq, tp, i)
        set_solution!(sol2, asol, i)
    end
    @test sol1.q[1:nt] == qs
    @test sol1.p[1:nt] == ps
    @test sol2.q[1:nt] == qs
    @test sol2.p[1:nt] == ps

    # sol1 = Solution(similar(pode, ics=(q=q1, p=p1)))
    # sol2 = Solution(similar(pode, ics=(q=q1, p=p1)))
    # for i = 1:nt
    #     for j = 1:ni
    #         tq .= Qs[i, j]
    #         tp .= Ps[i, j]
    #         asol.q .= Qs[i, j]
    #         asol.p .= Ps[i, j]
    #         set_solution!(sol1, tq, tp, i, j)
    #         set_solution!(sol2, asol, i, j)
    #     end
    # end
    # @test sol1.q[1:nt, :] == Qs
    # @test sol1.p[1:nt, :] == Ps
    # @test sol2.q[1:nt, :] == Qs
    # @test sol2.p[1:nt, :] == Ps

    # test nsave and nwrite parameters
    sol = Solution(similar(pode, tspan = 2 .* tspan(pode)), nsave = 2)
    @test sol.nt == 10

    sol = Solution(similar(pode, tspan = 2 .* tspan(pode)), nsave = 2, nwrite = 10)
    @test sol.nt == 5

    # test reset
    sol = Solution(pode)
    reset!(sol)
    @test sol.t[0] == t1
    @test sol.t[end] == t2
    @test offset(sol) == nt
end


@testset "$(rpad("DAE Solution",80))" begin
    asol = AtomicSolution(dae)

    # test constructors and general functionality
    sol = Solution(dae)
    @test typeof(sol) <: SolutionDAE

    # test_interface(sol) # TODO reactivate

    sol0 = Solution(similar(dae, ics=(q=z0, λ=λ0)))
    @test typeof(sol0) <: SolutionDAE

    # sol1 = Solution(similar(dae, ics=(q=z1, λ=λ1)))
    # @test typeof(sol1) <: SolutionDAE

    @test sol != sol0
    # @test sol != sol1

    @test nsave(sol) == 1
    @test ntime(sol) == nt
    @test timesteps(sol) == Δt .* collect(0:nt)

    # test initial conditions
    set_initial_conditions!(sol, t0, z0, λ0)
    get_initial_conditions!(sol, tz, tλ, 1)
    get_initial_conditions!(sol, asol, 1)

    @test tz == z0
    @test tλ == λ0
    @test asol.t == t0
    @test asol.q == z0
    @test asol.λ == λ0

    δt, δz, δλ = get_initial_conditions(sol, 1)
    @test δt == t0
    @test δz == z0
    @test δλ == λ0

    # set_initial_conditions!(sol1, similar(dae, t1, z2, λ2))
    # for j = 1:ni
    #     get_initial_conditions!(sol1, tz, tλ, j)
    #     @test tz == z2[j]
    #     @test tλ == λ2[j]
    # end

    # test set/get solution
    sol1 = Solution(similar(dae, ics=(q=z0, λ=λ0)))
    sol2 = Solution(similar(dae, ics=(q=z0, λ=λ0)))
    for i = 1:nt
        tz .= zs[i]
        tλ .= λs[i]
        asol.q .= zs[i]
        asol.λ .= λs[i]
        set_solution!(sol1, tz, tλ, i)
        set_solution!(sol2, asol, i)
    end
    @test sol1.q[1:nt] == zs
    @test sol1.λ[1:nt] == λs
    @test sol2.q[1:nt] == zs
    @test sol2.λ[1:nt] == λs

    # sol1 = Solution(similar(dae, ics=(q=z1, λ=λ1)))
    # sol2 = Solution(similar(dae, ics=(q=z1, λ=λ1)))
    # for i = 1:nt
    #     for j = 1:ni
    #         tz .= Zs[i, j]
    #         tλ .= Λs[i, j]
    #         asol.q .= Zs[i, j]
    #         asol.λ .= Λs[i, j]
    #         set_solution!(sol1, tz, tλ, i, j)
    #         set_solution!(sol2, asol, i, j)
    #     end
    # end
    # @test sol1.q[1:nt, :] == Zs
    # @test sol1.λ[1:nt, :] == Λs
    # @test sol2.q[1:nt, :] == Zs
    # @test sol2.λ[1:nt, :] == Λs

    # test nsave and nwrite parameters
    sol = Solution(similar(dae, tspan = 2 .* tspan(dae)), nsave = 2)
    @test sol.nt == 10

    sol = Solution(similar(dae, tspan = 2 .* tspan(dae)), nsave = 2, nwrite = 10)
    @test sol.nt == 5

    # test reset
    sol = Solution(dae)
    reset!(sol)
    @test sol.t[0] == t1
    @test sol.t[end] == t2
    @test offset(sol) == nt
end


@testset "$(rpad("PDAE Solution",80))" begin
    asol = AtomicSolution(pdae)

    # test constructors and general functionality
    sol = Solution(pdae)
    @test typeof(sol) <: SolutionPDAE

    # test_interface(sol) # TODO reactivate

    sol0 = Solution(similar(pdae, ics=(q=x0, p=y0, λ=μ0)))
    @test typeof(sol0) <: SolutionPDAE

    # sol1 = Solution(similar(pdae, ics=(q=x1, p=y1, λ=μ1)))
    # @test typeof(sol1) <: SolutionPDAE

    @test sol != sol0
    # @test sol != sol1

    @test nsave(sol) == 1
    @test ntime(sol) == nt
    @test timesteps(sol) == Δt .* collect(0:nt)

    # test initial conditions
    set_initial_conditions!(sol, t0, x0, y0, μ0)
    get_initial_conditions!(sol, tx, ty, tμ, 1)
    get_initial_conditions!(sol, asol, 1)

    @test tx == x0
    @test ty == y0
    @test tμ == μ0
    @test asol.t == t0
    @test asol.q == x0
    @test asol.p == y0
    @test asol.λ == μ0

    tx .= 0
    ty .= 0
    get_initial_conditions!(sol, tx, ty, 1)
    @test tx == x0
    @test ty == y0

    δt, δx, δy, δμ = get_initial_conditions(sol, 1)
    @test δt == t0
    @test δx == x0
    @test δy == y0
    @test δμ == μ0

    # set_initial_conditions!(sol1, similar(pdae, t1, x2, y2, μ2))
    # for j = 1:ni
    #     get_initial_conditions!(sol1, tx, ty, tμ, j)
    #     @test tx == x2[j]
    #     @test ty == y2[j]
    #     @test tμ == μ2[j]
    # end

    # test set/get solution
    sol1 = Solution(similar(pdae, ics=(q=x0, p=y0, λ=μ0)))
    sol2 = Solution(similar(pdae, ics=(q=x0, p=y0, λ=μ0)))
    sol3 = Solution(similar(pdae, ics=(q=x0, p=y0, λ=μ0)))
    for i = 1:nt
        tx .= xs[i]
        ty .= ys[i]
        tμ .= μs[i]
        asol.q .= xs[i]
        asol.p .= ys[i]
        asol.λ .= μs[i]
        set_solution!(sol1, tx, ty, tμ, i)
        set_solution!(sol2, tx, ty, i)
        set_solution!(sol3, asol, i)
    end
    @test sol1.q[1:nt] == xs
    @test sol1.p[1:nt] == ys
    @test sol1.λ[1:nt] == μs
    @test sol2.q[1:nt] == xs
    @test sol2.p[1:nt] == ys
    @test sol3.q[1:nt] == xs
    @test sol3.p[1:nt] == ys
    @test sol3.λ[1:nt] == μs

    # sol1 = Solution(similar(pdae, ics=(q=x1, p=y1, λ=μ1)))
    # sol2 = Solution(similar(pdae, ics=(q=x1, p=y1, λ=μ1)))
    # sol3 = Solution(similar(pdae, ics=(q=x1, p=y1, λ=μ1)))
    # for i = 1:nt
    #     for j = 1:ni
    #         tx .= Xs[i, j]
    #         ty .= Ys[i, j]
    #         tμ .= Ms[i, j]
    #         asol.q .= Xs[i, j]
    #         asol.p .= Ys[i, j]
    #         asol.λ .= Ms[i, j]
    #         set_solution!(sol1, tx, ty, tμ, i, j)
    #         set_solution!(sol2, tx, ty, i, j)
    #         set_solution!(sol3, asol, i, j)
    #     end
    # end
    # @test sol1.q[1:nt, :] == Xs
    # @test sol1.p[1:nt, :] == Ys
    # @test sol1.λ[1:nt, :] == Ms
    # @test sol2.q[1:nt, :] == Xs
    # @test sol2.p[1:nt, :] == Ys
    # @test sol3.q[1:nt, :] == Xs
    # @test sol3.p[1:nt, :] == Ys
    # @test sol3.λ[1:nt, :] == Ms

    # test nsave and nwrite parameters
    sol = Solution(similar(pdae, tspan = 2 .* tspan(pdae)), nsave = 2)
    @test sol.nt == 10

    sol = Solution(similar(pdae, tspan = 2 .* tspan(pdae)), nsave = 2, nwrite = 10)
    @test sol.nt == 5

    # test reset
    sol = Solution(pdae)
    reset!(sol)
    @test sol.t[0] == t1
    @test sol.t[end] == t2
    @test offset(sol) == nt
end
