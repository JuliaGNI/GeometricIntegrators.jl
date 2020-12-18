
using GeometricIntegrators.Config
using GeometricIntegrators.Integrators
using GeometricIntegrators.Integrators.Stochastic
using GeometricIntegrators.Solutions
using GeometricIntegrators.Tableaus
using GeometricProblems.KuboOscillator
using Test

using GeometricProblems.KuboOscillator: Δt, nt

include("utils.jl")

ode    = kubo_oscillator_ode()
sde1   = kubo_oscillator_sde_1()
sde2   = kubo_oscillator_sde_2()
sde3   = kubo_oscillator_sde_3()
psde1  = kubo_oscillator_psde_1()
psde2  = kubo_oscillator_psde_2()
psde3  = kubo_oscillator_psde_3()
spsde1 = kubo_oscillator_spsde_1()
spsde2 = kubo_oscillator_spsde_2()
spsde3 = kubo_oscillator_spsde_3()


@testset "$(rpad("Stochastic integrator constructors",80))" begin

    @test typeof(Integrator(sde1, TableauBurrageE1(), Δt)) <: IntegratorSERK
    @test typeof(Integrator(sde1, TableauStochasticGLRK(1), Δt)) <: IntegratorSIRK
    @test typeof(Integrator(sde1, TableauRosslerRS1(), Δt)) <: IntegratorWERK
    @test typeof(Integrator(sde1, TableauSRKw1(), Δt)) <: IntegratorWIRK
    @test typeof(Integrator(psde1, TableauStochasticStormerVerlet(), Δt)) <: IntegratorSIPRK
    @test typeof(Integrator(spsde1, TableauModifiedStochasticStormerVerlet(), Δt)) <: IntegratorSISPRK

end


@testset "$(rpad("SERK integrators",80))" begin

    int = Integrator(sde1, TableauBurrageE1(), Δt)
    sol = Solution(sde1, Δt, nt, conv=:strong)
    integrate!(int, sol)
    @test rel_energy_err(sol) < 2E-6

    int = Integrator(sde2, TableauBurrageE1(), Δt)
    sol = Solution(sde2, Δt, nt, conv=:strong)
    integrate!(int, sol)
    @test rel_energy_err(sol) < 2E-6

    int = Integrator(sde3, TableauBurrageE1(), Δt)
    sol = Solution(sde3, Δt, nt, conv=:strong)
    integrate!(int, sol)
    @test rel_energy_err(sol) < 2E-6

end


@testset "$(rpad("SIRK integrators",80))" begin

    int = Integrator(sde1, TableauStochasticGLRK(1), Δt)
    sol = Solution(sde1, Δt, nt, conv=:strong)
    integrate!(int, sol)
    @test rel_energy_err(sol) < 1E-14

    int = Integrator(sde2, TableauStochasticGLRK(1), Δt)
    sol = Solution(sde2, Δt, nt, conv=:strong)
    integrate!(int, sol)
    @test rel_energy_err(sol) < 1E-14

    int = Integrator(sde3, TableauStochasticGLRK(1), Δt)
    sol = Solution(sde3, Δt, nt, conv=:strong)
    integrate!(int, sol)
    @test rel_energy_err(sol) < 1E-14

    # check deterministic part of SIRK integrator to correspond to FIRK integrator
    int_sto = Integrator(sde1, TableauStochasticGLRK(1), Δt)
    sol_sto = Solution(sde1, Δt, zeros(1, nt), zeros(1, nt), nt, conv=:strong)
    integrate!(int_sto, sol_sto)

    int_det = Integrator(ode, TableauGLRK(1), Δt)
    sol_det = Solution(ode, Δt, nt)
    integrate!(int_det, sol_det)

    @test sol_sto.q[end] ≈ sol_det.q[end] atol=1E-14
end


@testset "$(rpad("WERK integrators",80))" begin

    int = Integrator(sde1, TableauRosslerRS1(), Δt)
    sol = Solution(sde1, Δt, nt, conv=:weak)
    integrate!(int, sol)
    @test rel_energy_err(sol) < 1E-5

    int = Integrator(sde2, TableauRosslerRS1(), Δt)
    sol = Solution(sde2, Δt, nt, conv=:weak)
    integrate!(int, sol)
    @test rel_energy_err(sol) < 1E-5

    int = Integrator(sde3, TableauRosslerRS1(), Δt)
    sol = Solution(sde3, Δt, nt, conv=:weak)
    integrate!(int, sol)
    @test rel_energy_err(sol) < 1E-5

end


@testset "$(rpad("WIRK integrators",80))" begin

    int = Integrator(sde1, TableauSRKw1(), Δt)
    sol = Solution(sde1, Δt, nt, conv=:weak)
    integrate!(int, sol)
    @test rel_energy_err(sol) < 1E-14

    int = Integrator(sde2, TableauSRKw1(), Δt)
    sol = Solution(sde2, Δt, nt, conv=:weak)
    integrate!(int, sol)
    @test rel_energy_err(sol) < 1E-14

    int = Integrator(sde3, TableauSRKw1(), Δt)
    sol = Solution(sde3, Δt, nt, conv=:weak)
    integrate!(int, sol)
    @test rel_energy_err(sol) < 1E-14

end


@testset "$(rpad("SIPRK integrators",80))" begin

    int = Integrator(psde1, TableauStochasticStormerVerlet(), Δt)
    sol = Solution(psde1, Δt, nt, conv=:strong)
    integrate!(int, sol)
    @test rel_energy_err(sol) < 1E-5

    int = Integrator(psde2, TableauStochasticStormerVerlet(), Δt)
    sol = Solution(psde2, Δt, nt, conv=:strong)
    integrate!(int, sol)
    @test rel_energy_err(sol) < 2E-5

    int = Integrator(psde3, TableauStochasticStormerVerlet(), Δt)
    sol = Solution(psde3, Δt, nt, conv=:strong)
    integrate!(int, sol)
    @test rel_energy_err(sol) < 2E-5

end


@testset "$(rpad("SISPRK integrators",80))" begin

    int = Integrator(spsde1, TableauModifiedStochasticStormerVerlet(), Δt)
    sol = Solution(spsde1, Δt, nt, conv=:strong)
    integrate!(int, sol)
    @test rel_energy_err(sol) < 2E-2

    int = Integrator(spsde2, TableauModifiedStochasticStormerVerlet(), Δt)
    sol = Solution(spsde2, Δt, nt, conv=:strong)
    integrate!(int, sol)
    @test rel_energy_err(sol) < 2E-2

    int = Integrator(spsde3, TableauModifiedStochasticStormerVerlet(), Δt)
    sol = Solution(spsde3, Δt, nt, conv=:strong)
    integrate!(int, sol)
    @test rel_energy_err(sol) < 2E-2

end
