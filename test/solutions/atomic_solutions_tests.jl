
using GeometricIntegrators.Common
using GeometricIntegrators.Solutions
using GeometricProblems.HarmonicOscillator
using Test


Δt = .1
t0 = 0.
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
    @test get_solution(asol) == (zero(t0), zero(x0))

    set_solution!(asol, (t0, x0))
    @test get_solution(asol) == (t0, x0)
    @test asol.t  == t0
    @test asol.q  == x0
    @test asol.t̄ == zero(t0)
    @test asol.q̄ == zero(x0)

    reset!(asol, Δt)
    @test asol.t̄ == t0
    @test asol.q̄ == x0

    update!(asol, v0)
    @test asol.t == t0  + Δt
    @test asol.q == x0 .+ v0

    set_solution!(asol, (t0, [2π,2π]))
    cut_periodic_solution!(asol, [2π, 0.])
    @test asol.q  == [0., 2π]

    set_solution!(asol, (t0, [-2π,2π]))
    cut_periodic_solution!(asol, [2π, 0.])
    @test asol.q  == [0., 2π]
end


@testset "$(rpad("Atomic PODE Solution",80))" begin
    asol = AtomicSolutionPODE(t0, q0, p0)
    @test get_solution(asol) == (zero(t0), zero(q0), zero(p0))

    set_solution!(asol, (t0, q0, p0))
    @test get_solution(asol) == (t0, q0, p0)
    @test asol.t  == t0
    @test asol.q  == q0
    @test asol.p  == p0
    @test asol.t̄ == zero(t0)
    @test asol.q̄ == zero(q0)
    @test asol.p̄ == zero(p0)

    reset!(asol, Δt)
    @test asol.t̄ == t0
    @test asol.q̄ == q0
    @test asol.p̄ == p0

    update!(asol, y0, z0)
    @test asol.t == t0  + Δt
    @test asol.q == q0 .+ y0
    @test asol.p == p0 .+ z0

    set_solution!(asol, (t0, [2π], [2π]))
    cut_periodic_solution!(asol, [2π])
    @test asol.q  == [0.]
    @test asol.p  == [2π]
end



@testset "$(rpad("Atomic DAE Solution",80))" begin
    asol = AtomicSolutionDAE(t0, x0, λ0)
    @test get_solution(asol) == (zero(t0), zero(x0), zero(λ0))

    set_solution!(asol, (t0, x0, λ0))
    @test get_solution(asol) == (t0, x0, λ0)
    @test asol.t  == t0
    @test asol.q  == x0
    @test asol.λ  == λ0
    @test asol.t̄ == zero(t0)
    @test asol.q̄ == zero(x0)
    @test asol.λ̄== zero(λ0)

    reset!(asol, Δt)
    @test asol.t̄ == t0
    @test asol.q̄ == x0
    @test asol.λ̄== λ0

    update!(asol, v0, λ1)
    @test asol.t == t0  + Δt
    @test asol.q == x0 .+ v0
    @test asol.λ == λ1
end



@testset "$(rpad("Atomic PDAE Solution",80))" begin
    asol = AtomicSolutionPDAE(t0, q0, p0, λ0)
    @test get_solution(asol) == (zero(t0), zero(q0), zero(p0), zero(λ0))

    set_solution!(asol, (t0, q0, p0, λ0))
    @test get_solution(asol) == (t0, q0, p0, λ0)
    @test asol.t  == t0
    @test asol.q  == q0
    @test asol.p  == p0
    @test asol.λ  == λ0
    @test asol.t̄ == zero(t0)
    @test asol.q̄ == zero(q0)
    @test asol.p̄ == zero(p0)
    @test asol.λ̄== zero(λ0)

    reset!(asol, Δt)
    @test asol.t̄ == t0
    @test asol.q̄ == q0
    @test asol.p̄ == p0
    @test asol.λ̄== λ0

    update!(asol, y0, z0, λ1)
    @test asol.t == t0  + Δt
    @test asol.q == q0 .+ y0
    @test asol.p == p0 .+ z0
    @test asol.λ == λ1
end
