
using GeometricBase
using GeometricIntegrators.Solutions
using GeometricProblems.HarmonicOscillator
using Test


Δt = .1
t0 = 0.
t1 = t0 + Δt
x0 = rand(2)
q0 = rand(1)
p0 = q0.^2
λ0 = rand(1)
λ1 = rand(1)
v0 = rand(2)
y0 = rand(1)
z0 = rand(1)
ΔW = rand(3)
ΔZ = rand(3)


@testset "$(rpad("Atomic Solution Constructors",80))" begin
    ode   = harmonic_oscillator_ode()
    dae   = harmonic_oscillator_dae()
    pode  = harmonic_oscillator_pode()
    pdae  = harmonic_oscillator_pdae()
    iode  = harmonic_oscillator_iode()
    idae  = harmonic_oscillator_idae()

    @test typeof(AtomicSolution(ode))   <: AtomicSolutionODE
    @test typeof(AtomicSolution(dae))   <: AtomicSolutionDAE
    @test typeof(AtomicSolution(pode))  <: AtomicSolutionPODE
    @test typeof(AtomicSolution(iode))  <: AtomicSolutionPODE
    @test typeof(AtomicSolution(pdae))  <: AtomicSolutionPDAE
    @test typeof(AtomicSolution(idae))  <: AtomicSolutionPDAE
end


@testset "$(rpad("Atomic ODE Solution",80))" begin
    asol = AtomicSolutionODE(t0, x0)

    @test asol.t  == t0
    @test asol.q  == x0
    @test asol.t̄ == zero(t0)
    @test asol.q̄ == zero(x0)

    @test current(asol) == (t = t0, q = x0)
    @test previous(asol) == (t = zero(t0), q = zero(x0))

    reset!(asol)
    @test current(asol) == (t = t0, q = x0)
    @test previous(asol) == (t = t0, q = x0)

    update!(asol, Δt, v0)
    @test asol.t == t0  + Δt
    @test asol.q == x0 .+ v0

    copy!(asol, (t = t1, q = [2π,2π]))
    @test asol.t == t1
    @test asol.q == [2π,2π]
    cut_periodic_solution!(asol, (q = [2π, 0.],))
    @test asol.q == [0.,2π]
end


@testset "$(rpad("Atomic PODE Solution",80))" begin
    asol = AtomicSolutionPODE(t0, q0, p0)

    @test asol.t  == t0
    @test asol.q  == q0
    @test asol.p  == p0
    @test asol.t̄ == zero(t0)
    @test asol.q̄ == zero(q0)
    @test asol.p̄ == zero(p0)

    @test current(asol) == (t = t0, q = q0, p = p0)
    @test previous(asol) == (t = zero(t0), q = zero(q0), p = zero(p0))

    reset!(asol)
    @test current(asol) == (t = t0, q = q0, p = p0)
    @test previous(asol) == (t = t0, q = q0, p = p0)

    update!(asol, Δt, y0, z0)
    @test asol.t == t0  + Δt
    @test asol.q == q0 .+ y0
    @test asol.p == p0 .+ z0

    copy!(asol, (t = t1, q = [2π], p = [2π]))
    @test asol.t == t1
    @test asol.q == [2π]
    @test asol.p == [2π]
    cut_periodic_solution!(asol, (q = [2π],))
    @test asol.q == [0.]
    @test asol.p == [2π]
end



@testset "$(rpad("Atomic DAE Solution",80))" begin
    asol = AtomicSolutionDAE(t0, x0, λ0)

    @test asol.t == t0
    @test asol.q == x0
    @test asol.λ == λ0
    @test asol.t̄ == zero(t0)
    @test asol.q̄ == zero(x0)
    @test asol.λ̄ == zero(λ0)

    @test current(asol) == (t = t0, q = x0, λ = λ0)
    @test previous(asol) == (t = zero(t0), q = zero(x0), λ = zero(λ0))

    reset!(asol)
    @test current(asol) == (t = t0, q = x0, λ = λ0)
    @test previous(asol) == (t = t0, q = x0, λ = λ0)

    update!(asol, Δt, v0, λ1)
    @test asol.t == t0  + Δt
    @test asol.q == x0 .+ v0
    @test asol.λ == λ1

    copy!(asol, (t = t1, q = [-2π,2π], λ = λ1))
    @test asol.t == t1
    @test asol.q == [-2π,2π]
    @test asol.λ == λ1
    cut_periodic_solution!(asol, (q = [2π, 0.],))
    @test asol.q == [0., 2π]
end



@testset "$(rpad("Atomic PDAE Solution",80))" begin
    asol = AtomicSolutionPDAE(t0, q0, p0, λ0)

    @test asol.t  == t0
    @test asol.q  == q0
    @test asol.p  == p0
    @test asol.λ  == λ0
    @test asol.t̄ == zero(t0)
    @test asol.q̄ == zero(q0)
    @test asol.p̄ == zero(p0)
    @test asol.λ̄ == zero(λ0)

    @test current(asol) == (t = t0, q = q0, p = p0, λ = λ0)
    @test previous(asol) == (t = zero(t0), q = zero(q0), p = zero(p0), λ = zero(λ0))

    reset!(asol)
    @test current(asol) == (t = t0, q = q0, p = p0, λ = λ0)
    @test previous(asol) == (t = t0, q = q0, p = p0, λ = λ0)

    update!(asol, Δt, y0, z0, λ1)
    @test asol.t == t0  + Δt
    @test asol.q == q0 .+ y0
    @test asol.p == p0 .+ z0
    @test asol.λ == λ1

    copy!(asol, (t = t1, q = [2π], p = [2π], λ = λ1))
    @test asol.t == t1
    @test asol.q == [2π]
    @test asol.p == [2π]
    @test asol.λ == λ1
    cut_periodic_solution!(asol, (q = [2π],))
    @test asol.q == [0.]
    @test asol.p == [2π]
end
