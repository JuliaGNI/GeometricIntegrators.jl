
using GeometricIntegrators.Config
using GeometricIntegrators.Integrators
using GeometricIntegrators.Tableaus
using GeometricIntegrators.Utils
using GeometricProblems.HarmonicOscillator
using Test

set_config(:nls_stol_break, 1E3)

using GeometricProblems.HarmonicOscillator: Δt, nt, refx

sode = harmonic_oscillator_sode()


sint = Integrator(sode, TableauLieA(), Δt)
ssol = integrate(sode, sint, nt)
@test rel_err(ssol.q, refx) < 5E-2

sintc = IntegratorComposition(sode, TableauLieA(), Δt)
ssolc = integrate(sode, sintc, nt)
@test ssol.q == ssolc.q

sint = Integrator(sode, TableauLieB(), Δt)
ssol = integrate(sode, sint, nt)
@test rel_err(ssol.q, refx) < 5E-2

sintc = IntegratorComposition(sode, TableauLieB(), Δt)
ssolc = integrate(sode, sintc, nt)
@test ssol.q == ssolc.q

sint = Integrator(sode, TableauStrang(), Δt)
ssol = integrate(sode, sint, nt)
@test rel_err(ssol.q, refx) < 1E-3

sintc = IntegratorComposition(sode, TableauStrang(), Δt)
ssolc = integrate(sode, sintc, nt)
@test ssol.q == ssolc.q

sint = Integrator(sode, TableauMcLachlan2(), Δt)
ssol = integrate(sode, sint, nt)
@test rel_err(ssol.q, refx) < 1E-4

sintc = IntegratorComposition(sode, TableauMcLachlan2(), Δt)
ssolc = integrate(sode, sintc, nt)
@test ssol.q == ssolc.q

sint = Integrator(sode, TableauMcLachlan4(), Δt)
ssol = integrate(sode, sint, nt)
@test rel_err(ssol.q, refx) < 5E-4

sintc = IntegratorComposition(sode, TableauMcLachlan4(), Δt)
ssolc = integrate(sode, sintc, nt)
@test ssol.q == ssolc.q

sint = Integrator(sode, TableauTripleJump(), Δt)
ssol = integrate(sode, sint, nt)
@test rel_err(ssol.q, refx) < 5E-6

sintc = IntegratorComposition(sode, TableauTripleJump(), Δt)
ssolc = integrate(sode, sintc, nt)
@test ssol.q == ssolc.q

sint = Integrator(sode, TableauSuzukiFractal(), Δt)
ssol = integrate(sode, sint, nt)
@test rel_err(ssol.q, refx) < 5E-7

sintc = IntegratorComposition(sode, TableauSuzukiFractal(), Δt)
ssolc = integrate(sode, sintc, nt)
@test ssol.q == ssolc.q


DT = eltype(sode.q₀)
D  = ndims(sode)

ints_glrk1 = (IntegratorConstructor(DT, D, TableauGLRK(1)), IntegratorConstructor(DT, D, TableauGLRK(1)))
ints_erk4  = (IntegratorConstructor(DT, D, TableauRK4()), IntegratorConstructor(DT, D, TableauRK4()))

sint = IntegratorComposition(sode, ints_erk4, TableauLieA(), Δt)
ssol = integrate(sode, sint, nt)
# println(rel_err(ssol.q, refx))
@test rel_err(ssol.q, refx) < 5E-2

sint = IntegratorComposition(sode, ints_glrk1, TableauStrang(), Δt)
ssol = integrate(sode, sint, nt)
# println(rel_err(ssol.q, refx))
@test rel_err(ssol.q, refx) < 1E-3
