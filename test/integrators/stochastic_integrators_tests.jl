
using GeometricIntegrators.Config
using GeometricIntegrators.Integrators
using GeometricIntegrators.Solutions
using GeometricIntegrators.Solvers
using GeometricIntegrators.Tableaus
using GeometricIntegrators.TestProblems.KuboOscillatorProblem
using Test

using GeometricIntegrators.TestProblems.KuboOscillatorProblem: Δt, nt

include("utils.jl")

sde1   = kubo_oscillator_sde_1()
sde2   = kubo_oscillator_sde_2()
sde3   = kubo_oscillator_sde_3()
psde1  = kubo_oscillator_psde_1()
psde2  = kubo_oscillator_psde_2()
psde3  = kubo_oscillator_psde_3()
spsde1  = kubo_oscillator_spsde_1()
spsde2  = kubo_oscillator_spsde_2()
spsde3  = kubo_oscillator_spsde_3()


@test typeof(Integrator(sde1, getTableauBurrageE1(), Δt)) <: IntegratorSERK
@test typeof(Integrator(sde1, getTableauStochasticGLRK(1), Δt)) <: IntegratorSIRK
@test typeof(Integrator(sde1, getTableauRosslerRS1(), Δt)) <: IntegratorWERK
@test typeof(Integrator(sde1, getTableauSRKw1(), Δt)) <: IntegratorWIRK
@test typeof(Integrator(psde1, getTableauStochasticStormerVerlet(), Δt)) <: IntegratorSIPRK
@test typeof(Integrator(spsde1, getTableauModifiedStochasticStormerVerlet(), Δt)) <: IntegratorSISPRK


@testset "$(rpad("SERK integrators",80))" begin
    int = Integrator(sde1, getTableauBurrageE1(), Δt)
    sol = Solution(sde1, Δt, nt, conv="strong")
    integrate!(int, sol)

    @test rel_energy_err_sde(sol) < 2E-6

    int = Integrator(sde2, getTableauBurrageE1(), Δt)
    sol = Solution(sde2, Δt, nt, conv="strong")
    integrate!(int, sol)

    @test rel_energy_err_sde(sol) < 2E-6

    int = Integrator(sde3, getTableauBurrageE1(), Δt)
    sol = Solution(sde3, Δt, nt, conv="strong")
    integrate!(int, sol)

    @test rel_energy_err_sde(sol) < 2E-6
end


@testset "$(rpad("SIRK integrators",80))" begin
    int = Integrator(sde1, getTableauStochasticGLRK(1), Δt)
    sol = Solution(sde1, Δt, nt, conv="strong")
    integrate!(int, sol)

    @test rel_energy_err_sde(sol) < 1E-14

    int = Integrator(sde2, getTableauStochasticGLRK(1), Δt)
    sol = Solution(sde2, Δt, nt, conv="strong")
    integrate!(int, sol)

    @test rel_energy_err_sde(sol) < 1E-14

    int = Integrator(sde3, getTableauStochasticGLRK(1), Δt)
    sol = Solution(sde3, Δt, nt, conv="strong")
    integrate!(int, sol)

    @test rel_energy_err_sde(sol) < 1E-14
end


@testset "$(rpad("WERK integrators",80))" begin
    int = Integrator(sde1, getTableauRosslerRS1(), Δt)
    sol = Solution(sde1, Δt, nt, conv="weak")
    integrate!(int, sol)

    @test rel_energy_err_sde(sol) < 1E-5

    int = Integrator(sde2, getTableauRosslerRS1(), Δt)
    sol = Solution(sde2, Δt, nt, conv="weak")
    integrate!(int, sol)

    @test rel_energy_err_sde(sol) < 1E-5

    int = Integrator(sde3, getTableauRosslerRS1(), Δt)
    sol = Solution(sde3, Δt, nt, conv="weak")
    integrate!(int, sol)

    @test rel_energy_err_sde(sol) < 1E-5
end


@testset "$(rpad("WIRK integrators",80))" begin
    int = Integrator(sde1, getTableauSRKw1(), Δt)
    sol = Solution(sde1, Δt, nt, conv="weak")
    integrate!(int, sol)

    @test rel_energy_err_sde(sol) < 1E-14

    int = Integrator(sde2, getTableauSRKw1(), Δt)
    sol = Solution(sde2, Δt, nt, conv="weak")
    integrate!(int, sol)

    @test rel_energy_err_sde(sol) < 1E-14

    int = Integrator(sde3, getTableauSRKw1(), Δt)
    sol = Solution(sde3, Δt, nt, conv="weak")
    integrate!(int, sol)

    @test rel_energy_err_sde(sol) < 1E-14
end


@testset "$(rpad("SIPRK integrators",80))" begin
    int = Integrator(psde1, getTableauStochasticStormerVerlet(), Δt)
    sol = Solution(psde1, Δt, nt, conv="strong")
    integrate!(int, sol)

    @test rel_energy_err_psde(sol) < 1E-5

    int = Integrator(psde2, getTableauStochasticStormerVerlet(), Δt)
    sol = Solution(psde2, Δt, nt, conv="strong")
    integrate!(int, sol)

    @test rel_energy_err_psde(sol) < 1E-5

    int = Integrator(psde3, getTableauStochasticStormerVerlet(), Δt)
    sol = Solution(psde3, Δt, nt, conv="strong")
    integrate!(int, sol)

    @test rel_energy_err_psde(sol) < 1E-5
end


@testset "$(rpad("SISPRK integrators",80))" begin
    int = Integrator(spsde1, getTableauModifiedStochasticStormerVerlet(), Δt)
    sol = Solution(spsde1, Δt, nt, conv="strong")
    integrate!(int, sol)

    @test rel_energy_err_psde(sol) < 0.02

    int = Integrator(spsde2, getTableauModifiedStochasticStormerVerlet(), Δt)
    sol = Solution(spsde2, Δt, nt, conv="strong")
    integrate!(int, sol)

    @test rel_energy_err_psde(sol) < 0.02

    int = Integrator(spsde3, getTableauModifiedStochasticStormerVerlet(), Δt)
    sol = Solution(spsde3, Δt, nt, conv="strong")
    integrate!(int, sol)

    @test rel_energy_err_psde(sol) < 0.02
end
