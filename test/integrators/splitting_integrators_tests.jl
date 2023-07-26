using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
using Test


sode = sodeproblem()
ref  = exact_solution(odeproblem())

sode1 = SubstepProblem(sode, one(timestep(sode)), 1)
sode2 = SubstepProblem(sode, one(timestep(sode)), 2)

@test_nowarn Integrator(sode1, ExactSolution())
@test_nowarn Integrator(sode2, ExactSolution())


ssol = integrate(sode, LieA())
@test relative_maximum_error(ssol, ref).q < 5E-2

ssolc = integrate(sode, Composition(LieA()))
@test ssol.q == ssolc.q

ssol = integrate(sode, LieB())
@test relative_maximum_error(ssol, ref).q < 5E-2

ssolc = integrate(sode, Composition(LieB()))
@test ssol.q == ssolc.q

ssol = integrate(sode, Strang())
@test relative_maximum_error(ssol, ref).q < 1E-3

ssolc = integrate(sode, Composition(Strang()))
@test ssol.q == ssolc.q

ssol = integrate(sode, StrangA())
@test relative_maximum_error(ssol, ref).q < 1E-3

ssolc = integrate(sode, Composition(StrangA()))
@test ssol.q == ssolc.q

ssol = integrate(sode, StrangB())
@test relative_maximum_error(ssol, ref).q < 1E-3

ssolc = integrate(sode, Composition(StrangB()))
@test ssol.q == ssolc.q


ssol1 = integrate(sode, Strang())
ssol2 = integrate(sode, StrangA())
ssol3 = integrate(sode, StrangB())
@test ssol1.q.d ≈ ssol2.q.d  atol=1e-14
@test ssol1.q.d ≈ ssol3.q.d  atol=1e-3


ssol = integrate(sode, McLachlan2())
@test relative_maximum_error(ssol, ref).q < 1E-4

ssolc = integrate(sode, Composition(McLachlan2()))
@test ssol.q == ssolc.q

ssol = integrate(sode, McLachlan4())
@test relative_maximum_error(ssol, ref).q < 5E-4

ssolc = integrate(sode, Composition(McLachlan4()))
@test ssol.q == ssolc.q

ssol = integrate(sode, TripleJump())
@test relative_maximum_error(ssol, ref).q < 5E-6

ssolc = integrate(sode, Composition(TripleJump()))
@test ssol.q == ssolc.q

ssol = integrate(sode, SuzukiFractal())
@test relative_maximum_error(ssol, ref).q < 5E-7

ssolc = integrate(sode, Composition(SuzukiFractal()))
@test ssol.q == ssolc.q


DT = datatype(sode)
D  = length(sode.ics.q)

ints_glrk1 = (Gauss(1), Gauss(1))
ints_erk4  = (RK4(), RK4())

sint = Composition(ints_erk4, LieA())
ssol = integrate(sode, sint)
# println(relative_maximum_error(ssol, ref).q)
@test relative_maximum_error(ssol, ref).q < 5E-2

sint = Composition(ints_glrk1, Strang())
ssol = integrate(sode, sint)
# println(relative_maximum_error(ssol, ref).q)
@test relative_maximum_error(ssol, ref).q < 1E-3
