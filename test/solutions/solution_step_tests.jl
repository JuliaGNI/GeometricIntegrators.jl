
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
    hode  = harmonic_oscillator_hode()
    hdae  = harmonic_oscillator_hdae()
    iode  = harmonic_oscillator_iode()
    idae  = harmonic_oscillator_idae()
    # lode  = harmonic_oscillator_lode()
    # ldae  = harmonic_oscillator_ldae()

    @test typeof(SolutionStep(ode))   <: SolutionStepODE
    @test typeof(SolutionStep(dae))   <: SolutionStepDAE
    @test typeof(SolutionStep(pode))  <: SolutionStepPODE
    @test typeof(SolutionStep(hode))  <: SolutionStepPODE
    @test typeof(SolutionStep(iode))  <: SolutionStepPODE
    # @test typeof(SolutionStep(lode))  <: SolutionStepPODE
    @test typeof(SolutionStep(pdae))  <: SolutionStepPDAE
    @test typeof(SolutionStep(hdae))  <: SolutionStepPDAE
    @test typeof(SolutionStep(idae))  <: SolutionStepPDAE
    # @test typeof(SolutionStep(ldae))  <: SolutionStepPDAE
end


@testset "$(rpad("Atomic ODE Solution",80))" begin
    solstep = SolutionStepODE(t0, x0)

    @test solstep.t  == t0
    @test solstep.q  == x0
    @test solstep.t̄ == zero(t0)
    @test solstep.q̄ == zero(x0)

    @test current(solstep) == (t = t0, q = x0)
    @test previous(solstep) == (t = zero(t0), q = zero(x0))

    reset!(solstep)
    @test current(solstep) == (t = t0, q = x0)
    @test previous(solstep) == (t = t0, q = x0)

    update!(solstep, Δt, v0)
    @test solstep.t == t0  + Δt
    @test solstep.q == x0 .+ v0

    copy!(solstep, (t = t1, q = [2π,2π]))
    @test solstep.t == t1
    @test solstep.q == [2π,2π]
    cut_periodic_solution!(solstep, (q = [2π, 0.],))
    @test solstep.q == [0.,2π]
end


@testset "$(rpad("Atomic PODE Solution",80))" begin
    solstep = SolutionStepPODE(t0, q0, p0)

    @test solstep.t  == t0
    @test solstep.q  == q0
    @test solstep.p  == p0
    @test solstep.t̄ == zero(t0)
    @test solstep.q̄ == zero(q0)
    @test solstep.p̄ == zero(p0)

    @test current(solstep) == (t = t0, q = q0, p = p0)
    @test previous(solstep) == (t = zero(t0), q = zero(q0), p = zero(p0))

    reset!(solstep)
    @test current(solstep) == (t = t0, q = q0, p = p0)
    @test previous(solstep) == (t = t0, q = q0, p = p0)

    update!(solstep, Δt, y0, z0)
    @test solstep.t == t0  + Δt
    @test solstep.q == q0 .+ y0
    @test solstep.p == p0 .+ z0

    copy!(solstep, (t = t1, q = [2π], p = [2π]))
    @test solstep.t == t1
    @test solstep.q == [2π]
    @test solstep.p == [2π]
    cut_periodic_solution!(solstep, (q = [2π],))
    @test solstep.q == [0.]
    @test solstep.p == [2π]
end



@testset "$(rpad("Atomic DAE Solution",80))" begin
    solstep = SolutionStepDAE(t0, x0, λ0)

    @test solstep.t == t0
    @test solstep.q == x0
    @test solstep.λ == λ0
    @test solstep.t̄ == zero(t0)
    @test solstep.q̄ == zero(x0)
    @test solstep.λ̄ == zero(λ0)

    @test current(solstep) == (t = t0, q = x0, λ = λ0)
    @test previous(solstep) == (t = zero(t0), q = zero(x0), λ = zero(λ0))

    reset!(solstep)
    @test current(solstep) == (t = t0, q = x0, λ = λ0)
    @test previous(solstep) == (t = t0, q = x0, λ = λ0)

    update!(solstep, Δt, v0, λ1)
    @test solstep.t == t0  + Δt
    @test solstep.q == x0 .+ v0
    @test solstep.λ == λ1

    copy!(solstep, (t = t1, q = [-2π,2π], λ = λ1))
    @test solstep.t == t1
    @test solstep.q == [-2π,2π]
    @test solstep.λ == λ1
    cut_periodic_solution!(solstep, (q = [2π, 0.],))
    @test solstep.q == [0., 2π]
end



@testset "$(rpad("Atomic PDAE Solution",80))" begin
    solstep = SolutionStepPDAE(t0, q0, p0, λ0)

    @test solstep.t  == t0
    @test solstep.q  == q0
    @test solstep.p  == p0
    @test solstep.λ  == λ0
    @test solstep.t̄ == zero(t0)
    @test solstep.q̄ == zero(q0)
    @test solstep.p̄ == zero(p0)
    @test solstep.λ̄ == zero(λ0)

    @test current(solstep) == (t = t0, q = q0, p = p0, λ = λ0)
    @test previous(solstep) == (t = zero(t0), q = zero(q0), p = zero(p0), λ = zero(λ0))

    reset!(solstep)
    @test current(solstep) == (t = t0, q = q0, p = p0, λ = λ0)
    @test previous(solstep) == (t = t0, q = q0, p = p0, λ = λ0)

    update!(solstep, Δt, y0, z0, λ1)
    @test solstep.t == t0  + Δt
    @test solstep.q == q0 .+ y0
    @test solstep.p == p0 .+ z0
    @test solstep.λ == λ1

    copy!(solstep, (t = t1, q = [2π], p = [2π], λ = λ1))
    @test solstep.t == t1
    @test solstep.q == [2π]
    @test solstep.p == [2π]
    @test solstep.λ == λ1
    cut_periodic_solution!(solstep, (q = [2π],))
    @test solstep.q == [0.]
    @test solstep.p == [2π]
end
