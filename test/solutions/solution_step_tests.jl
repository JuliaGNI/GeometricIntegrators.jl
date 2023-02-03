
using GeometricEquations.Tests.HarmonicOscillator
using GeometricIntegrators
using Test


Δt = .1
t0 = 0.
t1 = t0 + Δt
x0 = rand(2)
q0 = rand(1)
p0 = q0.^2
λ0 = rand(1)
λ1 = rand(1)
Δx = rand(2)
Δq = rand(1)
Δp = rand(1)
ΔW = rand(3)
ΔZ = rand(3)

ode   = odeproblem()
dae   = daeproblem()
pode  = podeproblem()
pdae  = pdaeproblem()
hode  = hodeproblem()
hdae  = hdaeproblem()
iode  = iodeproblem()
idae  = idaeproblem()
lode  = lodeproblem()
ldae  = ldaeproblem()


@testset "$(rpad("Solution Step Constructors",80))" begin
    @test typeof(SolutionStep(ode))   <: SolutionStepODE
    @test typeof(SolutionStep(dae))   <: SolutionStepDAE
    @test typeof(SolutionStep(pode))  <: SolutionStepPODE
    @test typeof(SolutionStep(hode))  <: SolutionStepPODE
    @test typeof(SolutionStep(iode))  <: SolutionStepPODE
    @test typeof(SolutionStep(lode))  <: SolutionStepPODE
    @test typeof(SolutionStep(pdae))  <: SolutionStepPDAE
    @test typeof(SolutionStep(hdae))  <: SolutionStepPDAE
    @test typeof(SolutionStep(idae))  <: SolutionStepPDAE
    @test typeof(SolutionStep(ldae))  <: SolutionStepPDAE
end


@testset "$(rpad("ODE Solution Step",80))" begin

    solstep = SolutionStepODE(t0, x0)

    @test solstep.t == solstep.t̄[0] == current(solstep).t
    @test solstep.q == solstep.q̄[0] == current(solstep).q

    @test solstep.t̄ == history(solstep).t
    @test solstep.q̄ == history(solstep).q

    solstep.t̄[0]  = t0
    solstep.q̄[0] .= x0

    solstep.t̄[1]  = zero(t0)
    solstep.q̄[1] .= zero(x0)

    @test current(solstep) == (t = t0, q = x0)
    @test previous(solstep) == (t = zero(t0), q = zero(x0))

    reset!(solstep, Δt)
    @test current(solstep) == (t = t0 + Δt, q = x0)
    @test previous(solstep) == (t = t0, q = x0)

    update!(solstep, Δx)
    @test solstep.t == t0  + Δt
    @test solstep.q == x0 .+ Δx

    copy!(solstep, (t = t1, q = [2π,2π]))
    @test solstep.t == t1
    @test solstep.q == [2π,2π]
    cut_periodic_solution!(solstep, (q = [2π, 0.],))
    @test solstep.q == [0.,2π]


    solstep = SolutionStep(ode)

    k = HarmonicOscillator.k
    ω = HarmonicOscillator.ω

    A = sqrt(solstep.q[2]^2 / k + solstep.q[1]^2)
    ϕ = asin(solstep.q[1] / A)

    @test solstep.t == initial_conditions(ode).t
    @test solstep.q == initial_conditions(ode).q

    @test solstep.t̄[1] ≈ - timestep(ode)
    @test solstep.q̄[1] ≈ [A * sin(- ω * timestep(ode) + ϕ), ω * A * cos(- ω * timestep(ode) + ϕ)]

end


@testset "$(rpad("PODE Solution Step",80))" begin

    solstep = SolutionStepPODE(t0, q0, p0)

    @test solstep.t == solstep.t̄[0] == current(solstep).t
    @test solstep.q == solstep.q̄[0] == current(solstep).q
    @test solstep.p == solstep.p̄[0] == current(solstep).p

    @test solstep.t̄ == history(solstep).t
    @test solstep.q̄ == history(solstep).q
    @test solstep.p̄ == history(solstep).p

    solstep.t̄[0]  = t0
    solstep.q̄[0] .= q0
    solstep.p̄[0] .= p0

    solstep.t̄[1]  = zero(t0)
    solstep.q̄[1] .= zero(q0)
    solstep.p̄[1] .= zero(p0)

    @test current(solstep) == (t = t0, q = q0, p = p0)
    @test previous(solstep) == (t = zero(t0), q = zero(q0), p = zero(p0))

    reset!(solstep, Δt)
    @test current(solstep) == (t = t0 + Δt, q = q0, p = p0)
    @test previous(solstep) == (t = t0, q = q0, p = p0)

    update!(solstep, Δq, Δp)
    @test solstep.t == t0  + Δt
    @test solstep.q == q0 .+ Δq
    @test solstep.p == p0 .+ Δp

    copy!(solstep, (t = t1, q = [2π], p = [2π]))
    @test solstep.t == t1
    @test solstep.q == [2π]
    @test solstep.p == [2π]
    cut_periodic_solution!(solstep, (q = [2π],))
    @test solstep.q == [0.]
    @test solstep.p == [2π]


    solstep = SolutionStep(pode)

    k = HarmonicOscillator.k
    ω = HarmonicOscillator.ω

    A = sqrt(solstep.p[1]^2 / k + solstep.q[1]^2)
    ϕ = asin(solstep.q[1] / A)

    @test solstep.t == initial_conditions(pode).t
    @test solstep.q == initial_conditions(pode).q
    @test solstep.p == initial_conditions(pode).p

    @test solstep.t̄[1] ≈ - timestep(pode)
    @test solstep.q̄[1] ≈ [A * sin(- ω * timestep(pode) + ϕ)]
    @test solstep.p̄[1] ≈ [ω * A * cos(- ω * timestep(pode) + ϕ)]

end


@testset "$(rpad("DAE Solution Step",80))" begin

    solstep = SolutionStepDAE(t0, x0, λ0)

    @test solstep.t == solstep.t̄[0] == current(solstep).t
    @test solstep.q == solstep.q̄[0] == current(solstep).q
    @test solstep.λ == solstep.λ̄[0] == current(solstep).λ

    @test solstep.t̄ == history(solstep).t
    @test solstep.q̄ == history(solstep).q
    @test solstep.λ̄ == history(solstep).λ

    solstep.t̄[0]  = t0
    solstep.q̄[0] .= x0
    solstep.λ̄[0] .= λ0

    solstep.t̄[1]  = zero(t0)
    solstep.q̄[1] .= zero(x0)
    solstep.λ̄[1] .= zero(λ0)

    @test current(solstep) == (t = t0, q = x0, λ = λ0)
    @test previous(solstep) == (t = zero(t0), q = zero(x0), λ = zero(λ0))

    reset!(solstep, Δt)
    @test current(solstep) == (t = t0 + Δt, q = x0, λ = λ0)
    @test previous(solstep) == (t = t0, q = x0, λ = λ0)

    update!(solstep, Δx, λ1)
    @test solstep.t == t0  + Δt
    @test solstep.q == x0 .+ Δx
    @test solstep.λ == λ1

    copy!(solstep, (t = t1, q = [-2π,2π], λ = λ1))
    @test solstep.t == t1
    @test solstep.q == [-2π,2π]
    @test solstep.λ == λ1
    cut_periodic_solution!(solstep, (q = [2π, 0.],))
    @test solstep.q == [0., 2π]


    solstep = SolutionStep(dae)

    k = HarmonicOscillator.k
    ω = HarmonicOscillator.ω

    A = sqrt(solstep.q[2]^2 / k + solstep.q[1]^2)
    ϕ = asin(solstep.q[1] / A)

    @test solstep.t == initial_conditions(dae).t
    @test solstep.q == initial_conditions(dae).q

    @test solstep.t̄[1] ≈ - timestep(dae)
    @test solstep.q̄[1][1:2] ≈ [A * sin(- ω * timestep(dae) + ϕ), ω * A * cos(- ω * timestep(dae) + ϕ)]

end


@testset "$(rpad("PDAE Solution Step",80))" begin

    solstep = SolutionStepPDAE(t0, q0, p0, λ0)

    @test solstep.t == solstep.t̄[0] == current(solstep).t
    @test solstep.q == solstep.q̄[0] == current(solstep).q
    @test solstep.p == solstep.p̄[0] == current(solstep).p
    @test solstep.λ == solstep.λ̄[0] == current(solstep).λ

    @test solstep.t̄ == history(solstep).t
    @test solstep.q̄ == history(solstep).q
    @test solstep.p̄ == history(solstep).p
    @test solstep.λ̄ == history(solstep).λ

    solstep.t̄[0]  = t0
    solstep.q̄[0] .= q0
    solstep.p̄[0] .= p0
    solstep.λ̄[0] .= λ0

    solstep.t̄[1]  = zero(t0)
    solstep.q̄[1] .= zero(q0)
    solstep.p̄[1] .= zero(p0)
    solstep.λ̄[1] .= zero(λ0)

    @test current(solstep) == (t = t0, q = q0, p = p0, λ = λ0)
    @test previous(solstep) == (t = zero(t0), q = zero(q0), p = zero(p0), λ = zero(λ0))

    reset!(solstep, Δt)
    @test current(solstep) == (t = t0 + Δt, q = q0, p = p0, λ = λ0)
    @test previous(solstep) == (t = t0, q = q0, p = p0, λ = λ0)

    update!(solstep, Δq, Δp, λ1)
    @test solstep.t == t0  + Δt
    @test solstep.q == q0 .+ Δq
    @test solstep.p == p0 .+ Δp
    @test solstep.λ == λ1

    copy!(solstep, (t = t1, q = [2π], p = [2π], λ = λ1))
    @test solstep.t == t1
    @test solstep.q == [2π]
    @test solstep.p == [2π]
    @test solstep.λ == λ1
    cut_periodic_solution!(solstep, (q = [2π],))
    @test solstep.q == [0.]
    @test solstep.p == [2π]


    solstep = SolutionStep(pdae)

    k = HarmonicOscillator.k
    ω = HarmonicOscillator.ω

    A = sqrt(solstep.p[1]^2 / k + solstep.q[1]^2)
    ϕ = asin(solstep.q[1] / A)

    @test solstep.t == initial_conditions(pdae).t
    @test solstep.q == initial_conditions(pdae).q
    @test solstep.p == initial_conditions(pdae).p

    @test solstep.t̄[1] ≈ - timestep(pdae)
    @test solstep.q̄[1][1] ≈ A * sin(- ω * timestep(pdae) + ϕ)
    @test solstep.p̄[1][1] ≈ ω * A * cos(- ω * timestep(pdae) + ϕ)

end
