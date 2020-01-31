
using GeometricIntegrators.CommonFunctions
using GeometricIntegrators.Solutions
using GeometricIntegrators.TestProblems.HarmonicOscillatorProblem
using GeometricIntegrators.TestProblems.KuboOscillatorProblem
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
    sde   = kubo_oscillator_sde_1()
    psde  = kubo_oscillator_psde_1()
    spsde = kubo_oscillator_spsde_1()

    @test typeof(AtomicSolution(ode))   <: AtomicSolutionODE
    @test typeof(AtomicSolution(dae))   <: AtomicSolutionDAE
    @test typeof(AtomicSolution(pode))  <: AtomicSolutionPODE
    @test typeof(AtomicSolution(iode))  <: AtomicSolutionPODE
    @test typeof(AtomicSolution(pdae))  <: AtomicSolutionPDAE
    @test typeof(AtomicSolution(idae))  <: AtomicSolutionPDAE
    @test typeof(AtomicSolution(sde))   <: AtomicSolutionSDE
    @test typeof(AtomicSolution(psde))  <: AtomicSolutionPSDE
    @test typeof(AtomicSolution(spsde)) <: AtomicSolutionPSDE
end


@testset "$(rpad("Atomic ODE Solution",80))" begin
    asol = AtomicSolutionODE(eltype(x0), typeof(t0), 2)
    @test get_solution(asol) == (zero(t0), zero(x0))

    set_solution!(asol, (t0, x0))
    @test get_solution(asol) == (t0, x0)
    @test asol.t  == t0
    @test asol.q  == x0
    @test asol.t̅ == zero(t0)
    @test asol.q̅ == zero(x0)

    reset!(asol, Δt)
    @test asol.t̅ == t0
    @test asol.q̅ == x0

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
    asol = AtomicSolutionPODE(eltype(q0), typeof(t0), 1)
    @test get_solution(asol) == (zero(t0), zero(q0), zero(p0))

    set_solution!(asol, (t0, q0, p0))
    @test get_solution(asol) == (t0, q0, p0)
    @test asol.t  == t0
    @test asol.q  == q0
    @test asol.p  == p0
    @test asol.t̅ == zero(t0)
    @test asol.q̅ == zero(q0)
    @test asol.p̅ == zero(p0)

    reset!(asol, Δt)
    @test asol.t̅ == t0
    @test asol.q̅ == q0
    @test asol.p̅ == p0

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
    asol = AtomicSolutionDAE(eltype(x0), typeof(t0), 2, 1)
    @test get_solution(asol) == (zero(t0), zero(x0), zero(λ0))

    set_solution!(asol, (t0, x0, λ0))
    @test get_solution(asol) == (t0, x0, λ0)
    @test asol.t  == t0
    @test asol.q  == x0
    @test asol.λ  == λ0
    @test asol.t̅ == zero(t0)
    @test asol.q̅ == zero(x0)
    @test asol.λ̅ == zero(λ0)

    reset!(asol, Δt)
    @test asol.t̅ == t0
    @test asol.q̅ == x0
    @test asol.λ̅ == λ0

    update!(asol, v0, λ1)
    @test asol.t == t0  + Δt
    @test asol.q == x0 .+ v0
    @test asol.λ == λ1
end



@testset "$(rpad("Atomic PDAE Solution",80))" begin
    asol = AtomicSolutionPDAE(eltype(q0), typeof(t0), 1, 1)
    @test get_solution(asol) == (zero(t0), zero(q0), zero(p0), zero(λ0))

    set_solution!(asol, (t0, q0, p0, λ0))
    @test get_solution(asol) == (t0, q0, p0, λ0)
    @test asol.t  == t0
    @test asol.q  == q0
    @test asol.p  == p0
    @test asol.λ  == λ0
    @test asol.t̅ == zero(t0)
    @test asol.q̅ == zero(q0)
    @test asol.p̅ == zero(p0)
    @test asol.λ̅ == zero(λ0)

    reset!(asol, Δt)
    @test asol.t̅ == t0
    @test asol.q̅ == q0
    @test asol.p̅ == p0
    @test asol.λ̅ == λ0

    update!(asol, y0, z0, λ1)
    @test asol.t == t0  + Δt
    @test asol.q == q0 .+ y0
    @test asol.p == p0 .+ z0
    @test asol.λ == λ1
end



@testset "$(rpad("Atomic SDE Solution",80))" begin
    asol = AtomicSolutionSDE(eltype(x0), typeof(t0), 2, 3)
    @test get_solution(asol) == (zero(t0), zero(x0))

    set_solution!(asol, (t0, x0))
    @test get_solution(asol) == (t0, x0)
    @test asol.t  == t0
    @test asol.q  == x0
    @test asol.t̅ == zero(t0)
    @test asol.q̅ == zero(x0)

    set_increments!(asol, (ΔW, ΔZ))
    δW = zero(ΔW)
    δZ = zero(ΔZ)
    get_increments!(asol, δW, δZ)
    @test get_increments(asol) == (ΔW, ΔZ)
    @test asol.ΔW == ΔW
    @test asol.ΔZ == ΔZ
    @test δW == ΔW
    @test δZ == ΔZ

    reset!(asol, Δt)
    @test asol.t̅ == t0
    @test asol.q̅ == x0

    update!(asol, v0)
    @test asol.t == t0  + Δt
    @test asol.q == x0 .+ v0
end



@testset "$(rpad("Atomic PSDE Solution",80))" begin
    asol = AtomicSolutionPSDE(eltype(q0), typeof(t0), 1, 3)
    @test get_solution(asol) == (zero(t0), zero(q0), zero(p0))
    @test get_increments(asol) == (zero(ΔW), zero(ΔZ))

    set_solution!(asol, (t0, q0, p0))
    @test get_solution(asol) == (t0, q0, p0)
    @test asol.t  == t0
    @test asol.q  == q0
    @test asol.p  == p0
    @test asol.t̅  == zero(t0)
    @test asol.q̅  == zero(q0)
    @test asol.p̅  == zero(p0)

    set_increments!(asol, (ΔW, ΔZ))
    δW = zero(ΔW)
    δZ = zero(ΔZ)
    get_increments!(asol, δW, δZ)
    @test get_increments(asol) == (ΔW, ΔZ)
    @test asol.ΔW == ΔW
    @test asol.ΔZ == ΔZ
    @test δW == ΔW
    @test δZ == ΔZ

    reset!(asol, Δt)
    @test asol.t̅ == t0
    @test asol.q̅ == q0
    @test asol.p̅ == p0

    update!(asol, y0, z0)
    @test asol.t == t0  + Δt
    @test asol.q == q0 .+ y0
    @test asol.p == p0 .+ z0
end
