
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

sint = Integrator(sode, getTableauLieB(), Δt)
ssol = integrate(sode, sint, nt)
@test rel_err(ssol.q, refx) < 5E-2

sint = Integrator(sode, getTableauStrang(), Δt)
ssol = integrate(sode, sint, nt)
@test rel_err(ssol.q, refx) < 1E-3

sint = Integrator(sode, getTableauMcLachlan2(), Δt)
ssol = integrate(sode, sint, nt)
@test rel_err(ssol.q, refx) < 1E-4

sint = Integrator(sode, getTableauMcLachlan4(), Δt)
ssol = integrate(sode, sint, nt)
@test rel_err(ssol.q, refx) < 5E-4

sint = Integrator(sode, getTableauTripleJump(), Δt)
ssol = integrate(sode, sint, nt)
@test rel_err(ssol.q, refx) < 5E-6

sint = Integrator(sode, getTableauSuzukiFractal(), Δt)
ssol = integrate(sode, sint, nt)
@test rel_err(ssol.q, refx) < 5E-7
