
using GeometricIntegrators.Integrators
using GeometricIntegrators.Tableaus
using GeometricIntegrators.TestProblems.HarmonicOscillatorProblem
using GeometricIntegrators.Utils
using Test

using GeometricIntegrators.TestProblems.HarmonicOscillatorProblem: Δt, nt, refx

sode = harmonic_oscillator_sode()


sint = Integrator(sode, getTableauLieA(), Δt)
ssol = integrate(sode, sint, nt)
@test rel_err(ssol.q, refx) < 5E-2

sintc = IntegratorComposition(sode, getTableauLieA(), Δt)
ssolc = integrate(sode, sintc, nt)
@test ssol.q == ssolc.q

sint = Integrator(sode, getTableauLieB(), Δt)
ssol = integrate(sode, sint, nt)
@test rel_err(ssol.q, refx) < 5E-2

sintc = IntegratorComposition(sode, getTableauLieB(), Δt)
ssolc = integrate(sode, sintc, nt)
@test ssol.q == ssolc.q

sint = Integrator(sode, getTableauStrang(), Δt)
ssol = integrate(sode, sint, nt)
@test rel_err(ssol.q, refx) < 1E-3

sintc = IntegratorComposition(sode, getTableauStrang(), Δt)
ssolc = integrate(sode, sintc, nt)
@test ssol.q == ssolc.q

sint = Integrator(sode, getTableauMcLachlan2(), Δt)
ssol = integrate(sode, sint, nt)
@test rel_err(ssol.q, refx) < 1E-4

sintc = IntegratorComposition(sode, getTableauMcLachlan2(), Δt)
ssolc = integrate(sode, sintc, nt)
@test ssol.q == ssolc.q

sint = Integrator(sode, getTableauMcLachlan4(), Δt)
ssol = integrate(sode, sint, nt)
@test rel_err(ssol.q, refx) < 5E-4

sintc = IntegratorComposition(sode, getTableauMcLachlan4(), Δt)
ssolc = integrate(sode, sintc, nt)
@test ssol.q == ssolc.q

sint = Integrator(sode, getTableauTripleJump(), Δt)
ssol = integrate(sode, sint, nt)
@test rel_err(ssol.q, refx) < 5E-6

sintc = IntegratorComposition(sode, getTableauTripleJump(), Δt)
ssolc = integrate(sode, sintc, nt)
@test ssol.q == ssolc.q

sint = Integrator(sode, getTableauSuzukiFractal(), Δt)
ssol = integrate(sode, sint, nt)
@test rel_err(ssol.q, refx) < 5E-7

sintc = IntegratorComposition(sode, getTableauSuzukiFractal(), Δt)
ssolc = integrate(sode, sintc, nt)
@test ssol.q == ssolc.q
