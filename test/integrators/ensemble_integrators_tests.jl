using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
using Test


odes  = odeensemble()
podes = podeensemble()

refs  = exact_solution(odes)
prefs = exact_solution(podes)


@testset "$(rpad("ODE Ensembles",80))" begin

    sols = integrate(odes, Gauss(2))
    @test relative_maximum_error(sols, refs).q < 8E-6

end


@testset "$(rpad("PODE Ensembles",80))" begin

    psols = integrate(podes, Gauss(2))
    @test relative_maximum_error(psols, prefs).q < 6E-6
    @test relative_maximum_error(psols, prefs).p < 8E-6

end
