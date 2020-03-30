
using GeometricIntegrators.CommonFunctions
using GeometricIntegrators.Equations
using Test


@testset "$(rpad("Harmonic Oscillator",80))" begin

    using GeometricIntegrators.TestProblems.HarmonicOscillatorProblem

    ode  = harmonic_oscillator_ode()
    pode = harmonic_oscillator_pode()
    iode = harmonic_oscillator_iode()
    dae  = harmonic_oscillator_dae()
    idae = harmonic_oscillator_idae()
    pdae = harmonic_oscillator_pdae()

    ode_equs = get_function_tuple(ode)
    @test_nowarn ode_equs[:v](ode.t₀, ode.q₀, zero(ode.q₀))

    pode_equs = get_function_tuple(pode)
    @test_nowarn pode_equs[:v](pode.t₀, pode.q₀, pode.p₀, zero(pode.q₀))
    @test_nowarn pode_equs[:f](pode.t₀, pode.q₀, pode.p₀, zero(pode.p₀))

    iode_equs = get_function_tuple(iode)
    @test_nowarn iode_equs[:ϑ](iode.t₀, iode.q₀, iode.p₀, zero(iode.q₀))
    @test_nowarn iode_equs[:f](iode.t₀, iode.q₀, iode.p₀, zero(iode.p₀))
    @test_nowarn iode_equs[:v](iode.t₀, iode.q₀, zero(iode.q₀))

    dae_equs = get_function_tuple(dae)
    @test_nowarn dae_equs[:v](dae.t₀, dae.q₀, zero(dae.q₀))
    @test_nowarn dae_equs[:u](dae.t₀, dae.q₀, dae.λ₀, zero(dae.q₀))
    @test_nowarn dae_equs[:ϕ](dae.t₀, dae.q₀, zero(dae.λ₀))

    pdae_equs = get_function_tuple(pdae)
    @test_nowarn pdae_equs[:v](pdae.t₀, pdae.q₀, pdae.p₀, zero(pdae.q₀))
    @test_nowarn pdae_equs[:f](pdae.t₀, pdae.q₀, pdae.p₀, zero(pdae.p₀))
    @test_nowarn pdae_equs[:u](pdae.t₀, pdae.q₀, pdae.p₀, pdae.λ₀, zero(pdae.q₀))
    @test_nowarn pdae_equs[:g](pdae.t₀, pdae.q₀, pdae.p₀, pdae.λ₀, zero(pdae.p₀))
    @test_nowarn pdae_equs[:ϕ](pdae.t₀, pdae.q₀, pdae.p₀, zero(pdae.λ₀))

    idae_equs = get_function_tuple(idae)
    @test_nowarn idae_equs[:ϑ](idae.t₀, idae.q₀, idae.λ₀, zero(idae.q₀))
    @test_nowarn idae_equs[:f](idae.t₀, idae.q₀, idae.λ₀, zero(idae.p₀))
    @test_nowarn idae_equs[:u](idae.t₀, idae.q₀, pdae.p₀, idae.λ₀, zero(idae.q₀))
    @test_nowarn idae_equs[:g](idae.t₀, idae.q₀, pdae.p₀, idae.λ₀, zero(idae.p₀))
    @test_nowarn idae_equs[:ϕ](idae.t₀, idae.q₀, idae.p₀, zero(idae.λ₀))
    @test_nowarn idae_equs[:v](idae.t₀, idae.q₀, zero(idae.q₀))

end


@testset "$(rpad("Lotka-Volterra 2d",80))" begin

    using GeometricIntegrators.TestProblems.LotkaVolterra2dProblem

    ode  = lotka_volterra_2d_ode()
    pode = lotka_volterra_2d_pode()
    iode = lotka_volterra_2d_iode()
    vode = lotka_volterra_2d_vode()
    pdae = lotka_volterra_2d_pdae()
    idae = lotka_volterra_2d_idae()
    vdae = lotka_volterra_2d_vdae()

    ode_equs = get_function_tuple(ode)
    @test_nowarn ode_equs[:v](ode.t₀, ode.q₀, zero(ode.q₀))
    @test_nowarn ode_equs[:h](ode.t₀, ode.q₀)

    pode_equs = get_function_tuple(pode)
    @test_nowarn pode_equs[:v](pode.t₀, pode.q₀, pode.p₀, zero(pode.q₀))
    @test_nowarn pode_equs[:f](pode.t₀, pode.q₀, pode.p₀, zero(pode.p₀))
    @test_nowarn pode_equs[:h](pode.t₀, pode.q₀, pode.p₀)

    iode_equs = get_function_tuple(iode)
    @test_nowarn iode_equs[:ϑ](iode.t₀, iode.q₀, iode.p₀, zero(iode.q₀))
    @test_nowarn iode_equs[:f](iode.t₀, iode.q₀, iode.p₀, zero(iode.p₀))
    @test_nowarn iode_equs[:v](iode.t₀, iode.q₀, zero(iode.q₀))
    @test_nowarn iode_equs[:h](iode.t₀, iode.q₀)

    vode_equs = get_function_tuple(vode)
    @test_nowarn vode_equs[:ϑ](vode.t₀, vode.q₀, vode.p₀, zero(vode.q₀))
    @test_nowarn vode_equs[:f](vode.t₀, vode.q₀, vode.p₀, zero(vode.p₀))
    @test_nowarn vode_equs[:v](vode.t₀, vode.q₀, zero(vode.q₀))
    @test_nowarn vode_equs[:h](vode.t₀, vode.q₀)

    pdae_equs = get_function_tuple(pdae)
    @test_nowarn pdae_equs[:v](pdae.t₀, pdae.q₀, pdae.p₀, zero(pdae.q₀))
    @test_nowarn pdae_equs[:f](pdae.t₀, pdae.q₀, pdae.p₀, zero(pdae.p₀))
    @test_nowarn pdae_equs[:u](pdae.t₀, pdae.q₀, pdae.p₀, pdae.λ₀, zero(pdae.q₀))
    @test_nowarn pdae_equs[:g](pdae.t₀, pdae.q₀, pdae.p₀, pdae.λ₀, zero(pdae.p₀))
    @test_nowarn pdae_equs[:ϕ](pdae.t₀, pdae.q₀, pdae.p₀, zero(pdae.λ₀))
    @test_nowarn pdae_equs[:h](pdae.t₀, pdae.q₀, pdae.p₀)

    idae_equs = get_function_tuple(idae)
    @test_nowarn idae_equs[:ϑ](idae.t₀, idae.q₀, idae.p₀, zero(idae.q₀))
    @test_nowarn idae_equs[:f](idae.t₀, idae.q₀, idae.p₀, zero(idae.p₀))
    @test_nowarn idae_equs[:u](idae.t₀, idae.q₀, idae.p₀, idae.λ₀, zero(idae.q₀))
    @test_nowarn idae_equs[:g](idae.t₀, idae.q₀, idae.p₀, idae.λ₀, zero(idae.p₀))
    @test_nowarn idae_equs[:ϕ](idae.t₀, idae.q₀, idae.p₀, zero(idae.λ₀))
    @test_nowarn idae_equs[:v](idae.t₀, idae.q₀, zero(idae.q₀))
    @test_nowarn idae_equs[:h](idae.t₀, idae.q₀)

    vdae_equs = get_function_tuple(vdae)
    @test_nowarn vdae_equs[:ϑ](vdae.t₀, vdae.q₀, vdae.λ₀, zero(vdae.q₀))
    @test_nowarn vdae_equs[:f](vdae.t₀, vdae.q₀, vdae.λ₀, zero(vdae.p₀))
    @test_nowarn vdae_equs[:g](vdae.t₀, vdae.q₀, vdae.λ₀, zero(vdae.p₀))
    @test_nowarn vdae_equs[:g̅](vdae.t₀, vdae.q₀, vdae.λ₀, zero(vdae.p₀))
    @test_nowarn vdae_equs[:ϕ](vdae.t₀, vdae.q₀, vdae.p₀, zero(vdae.λ₀))
    @test_nowarn vdae_equs[:ψ](vdae.t₀, vdae.q₀, vdae.p₀, zero(vdae.q₀), zero(vdae.p₀), zero(vdae.λ₀))
    @test_nowarn vdae_equs[:v](vdae.t₀, vdae.q₀, zero(vdae.q₀))
    @test_nowarn vdae_equs[:h](vdae.t₀, vdae.q₀)

end


@testset "$(rpad("Kubo Oscillator",80))" begin

    using GeometricIntegrators.TestProblems.KuboOscillatorProblem

    sde   = kubo_oscillator_sde_1()
    psde  = kubo_oscillator_psde_1()
    spsde = kubo_oscillator_spsde_1()

    sde_equs = get_function_tuple(sde)
    @test_nowarn sde_equs[:v](sde.t₀, sde.q₀, zero(sde.q₀))
    @test_nowarn sde_equs[:B](sde.t₀, sde.q₀, zeros(eltype(sde.q₀), sde.d, sde.m))
    @test_nowarn sde_equs[:B](sde.t₀, sde.q₀, zeros(eltype(sde.q₀), sde.d, sde.m), 1)

    psde_equs = get_function_tuple(psde)
    @test_nowarn psde_equs[:v](psde.t₀, psde.q₀, psde.p₀, zero(psde.q₀))
    @test_nowarn psde_equs[:f](psde.t₀, psde.q₀, psde.p₀, zero(psde.p₀))
    @test_nowarn psde_equs[:B](psde.t₀, psde.q₀, psde.p₀, zero(psde.q₀))
    @test_nowarn psde_equs[:G](psde.t₀, psde.q₀, psde.p₀, zero(psde.p₀))
    @test_nowarn psde_equs[:B](psde.t₀, psde.q₀, psde.p₀, zeros(eltype(psde.q₀), psde.d, psde.m))
    @test_nowarn psde_equs[:G](psde.t₀, psde.q₀, psde.p₀, zeros(eltype(psde.p₀), psde.d, psde.m))

    spsde_equs = get_function_tuple(spsde)
    @test_nowarn spsde_equs[:v](spsde.t₀, spsde.q₀, spsde.p₀, zero(spsde.q₀))
    @test_nowarn spsde_equs[:f1](spsde.t₀, spsde.q₀, spsde.p₀, zero(spsde.p₀))
    @test_nowarn spsde_equs[:f2](spsde.t₀, spsde.q₀, spsde.p₀, zero(spsde.p₀))
    @test_nowarn spsde_equs[:B](spsde.t₀, spsde.q₀, spsde.p₀, zero(spsde.q₀))
    @test_nowarn spsde_equs[:G1](spsde.t₀, spsde.q₀, spsde.p₀, zero(spsde.p₀))
    @test_nowarn spsde_equs[:G2](spsde.t₀, spsde.q₀, spsde.p₀, zero(spsde.p₀))
    @test_nowarn spsde_equs[:B](spsde.t₀, spsde.q₀, spsde.p₀, zeros(eltype(spsde.q₀), spsde.d, spsde.m))
    @test_nowarn spsde_equs[:G1](spsde.t₀, spsde.q₀, spsde.p₀, zeros(eltype(spsde.p₀), spsde.d, spsde.m))
    @test_nowarn spsde_equs[:G2](spsde.t₀, spsde.q₀, spsde.p₀, zeros(eltype(spsde.p₀), spsde.d, spsde.m))

end
