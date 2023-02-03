using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
using SimpleSolvers
using Test

using GeometricProblems.HarmonicOscillator: reference_solution

sode = harmonic_oscillator_sode()


ssol = integrate(sode, LieA())
@test relative_maximum_error(ssol.q, reference_solution) < 5E-2

# ssolc = integrate(sode, Composition(sode, LieA()))
# @test ssol.q == ssolc.q

ssol = integrate(sode, LieB())
@test relative_maximum_error(ssol.q, reference_solution) < 5E-2

# ssolc = integrate(sode, Composition(sode, LieB()))
# @test ssol.q == ssolc.q

ssol = integrate(sode, Strang())
@test relative_maximum_error(ssol.q, reference_solution) < 1E-3

# ssolc = integrate(sode, Composition(sode, Strang()))
# @test ssol.q == ssolc.q

ssol = integrate(sode, StrangA())
@test relative_maximum_error(ssol.q, reference_solution) < 1E-3

# ssolc = integrate(sode, Composition(sode, StrangA()))
# @test ssol.q == ssolc.q

ssol = integrate(sode, StrangB())
@test relative_maximum_error(ssol.q, reference_solution) < 1E-3

# ssolc = integrate(sode, Composition(sode, StrangB()))
# @test ssol.q == ssolc.q

ssol1 = integrate(sode, Strang())
ssol2 = integrate(sode, StrangA())
ssol3 = integrate(sode, StrangB())
@test ssol1.q.d ≈ ssol2.q.d  atol=1e-14
@test ssol1.q.d ≈ ssol3.q.d  atol=1e-3

ssol = integrate(sode, McLachlan2())
@test relative_maximum_error(ssol.q, reference_solution) < 1E-4

# ssolc = integrate(sode, Composition(sode, McLachlan2()))
# @test ssol.q == ssolc.q

ssol = integrate(sode, McLachlan4())
@test relative_maximum_error(ssol.q, reference_solution) < 5E-4

# ssolc = integrate(sode, Composition(sode, McLachlan4()))
# @test ssol.q == ssolc.q

ssol = integrate(sode, TripleJump())
@test relative_maximum_error(ssol.q, reference_solution) < 5E-6

# ssolc = integrate(sode, Composition(sode, TripleJump()))
# @test ssol.q == ssolc.q

ssol = integrate(sode, SuzukiFractal())
@test relative_maximum_error(ssol.q, reference_solution) < 5E-7

# ssolc = integrate(sode, Composition(sode, SuzukiFractal()))
# @test ssol.q == ssolc.q


# DT = datatype(sode)
# D  = length(sode.ics.q)

# ints_glrk1 = (IntegratorConstructor(DT, D, TableauGauss(1)), IntegratorConstructor(DT, D, TableauGauss(1)))
# ints_erk4  = (IntegratorConstructor(DT, D, TableauRK4()), IntegratorConstructor(DT, D, TableauRK4()))

# sint = IntegratorComposition(sode, ints_erk4, LieA())
# ssol = integrate(sode, sint)
# # println(relative_maximum_error(ssol.q, reference_solution))
# @test relative_maximum_error(ssol.q, reference_solution) < 5E-2

# sint = IntegratorComposition(sode, ints_glrk1, Strang())
# ssol = integrate(sode, sint)
# # println(relative_maximum_error(ssol.q, reference_solution))
# @test relative_maximum_error(ssol.q, reference_solution) < 1E-3
