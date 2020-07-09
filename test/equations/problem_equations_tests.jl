
using GeometricIntegrators.Common
using GeometricIntegrators.Equations
using Test


@testset "$(rpad("Harmonic Oscillator",80))" begin

    using GeometricProblems.HarmonicOscillator

    ode  = harmonic_oscillator_ode()
    pode = harmonic_oscillator_pode()
    iode = harmonic_oscillator_iode()
    dae  = harmonic_oscillator_dae()
    idae = harmonic_oscillator_idae()
    pdae = harmonic_oscillator_pdae()

    ode_equs = get_function_tuple(ode)
    @test_nowarn ode_equs[:v](ode.t₀, ode.q₀[begin], zero(ode.q₀[begin]))

    iode_equs = get_function_tuple(iode)
    @test_nowarn iode_equs[:ϑ](iode.t₀, iode.q₀[begin], iode.p₀[begin], zero(iode.q₀[begin]))
    @test_nowarn iode_equs[:f](iode.t₀, iode.q₀[begin], iode.p₀[begin], zero(iode.p₀[begin]))
    @test_nowarn iode_equs[:v̄](iode.t₀, iode.q₀[begin], zero(iode.q₀[begin]))
    @test_nowarn iode_equs[:f̄](iode.t₀, iode.q₀[begin], iode.p₀[begin], zero(iode.q₀[begin]))

    pode_equs = get_function_tuple(pode)
    @test_nowarn pode_equs[:v](pode.t₀, pode.q₀[begin], pode.p₀[begin], zero(pode.q₀[begin]))
    @test_nowarn pode_equs[:f](pode.t₀, pode.q₀[begin], pode.p₀[begin], zero(pode.p₀[begin]))

    dae_equs = get_function_tuple(dae)
    @test_nowarn dae_equs[:v](dae.t₀, dae.q₀[begin], zero(dae.q₀[begin]))
    @test_nowarn dae_equs[:u](dae.t₀, dae.q₀[begin], dae.λ₀[begin], zero(dae.q₀[begin]))
    @test_nowarn dae_equs[:ϕ](dae.t₀, dae.q₀[begin], zero(dae.λ₀[begin]))
    @test_nowarn dae_equs[:v̄](dae.t₀, dae.q₀[begin], zero(dae.q₀[begin]))

    idae_equs = get_function_tuple(idae)
    @test_nowarn idae_equs[:ϑ](idae.t₀, idae.q₀[begin], idae.λ₀[begin], zero(idae.q₀[begin]))
    @test_nowarn idae_equs[:f](idae.t₀, idae.q₀[begin], idae.λ₀[begin], zero(idae.p₀[begin]))
    @test_nowarn idae_equs[:u](idae.t₀, idae.q₀[begin], pdae.p₀[begin], idae.λ₀[begin], zero(idae.q₀[begin]))
    @test_nowarn idae_equs[:g](idae.t₀, idae.q₀[begin], pdae.p₀[begin], idae.λ₀[begin], zero(idae.p₀[begin]))
    @test_nowarn idae_equs[:ϕ](idae.t₀, idae.q₀[begin], idae.p₀[begin], zero(idae.λ₀[begin]))
    @test_nowarn idae_equs[:v̄](idae.t₀, idae.q₀[begin], zero(idae.q₀[begin]))
    @test_nowarn idae_equs[:f̄](idae.t₀, idae.q₀[begin], idae.p₀[begin], zero(idae.q₀[begin]))

    pdae_equs = get_function_tuple(pdae)
    @test_nowarn pdae_equs[:v](pdae.t₀, pdae.q₀[begin], pdae.p₀[begin], zero(pdae.q₀[begin]))
    @test_nowarn pdae_equs[:f](pdae.t₀, pdae.q₀[begin], pdae.p₀[begin], zero(pdae.p₀[begin]))
    @test_nowarn pdae_equs[:u](pdae.t₀, pdae.q₀[begin], pdae.p₀[begin], pdae.λ₀[begin], zero(pdae.q₀[begin]))
    @test_nowarn pdae_equs[:g](pdae.t₀, pdae.q₀[begin], pdae.p₀[begin], pdae.λ₀[begin], zero(pdae.p₀[begin]))
    @test_nowarn pdae_equs[:ϕ](pdae.t₀, pdae.q₀[begin], pdae.p₀[begin], zero(pdae.λ₀[begin]))
    @test_nowarn pdae_equs[:v̄](pdae.t₀, pdae.q₀[begin], pdae.p₀[begin], zero(pdae.q₀[begin]))
    @test_nowarn pdae_equs[:f̄](pdae.t₀, pdae.q₀[begin], pdae.p₀[begin], zero(pdae.p₀[begin]))

end


@testset "$(rpad("Lotka-Volterra 2d",80))" begin

    using GeometricProblems.LotkaVolterra2d

    ode  = lotka_volterra_2d_ode()
    hode = lotka_volterra_2d_hode()
    iode = lotka_volterra_2d_iode()
    pode = lotka_volterra_2d_pode()
    vode = lotka_volterra_2d_vode()
    dae  = lotka_volterra_2d_dae()
    hdae = lotka_volterra_2d_hdae()
    idae = lotka_volterra_2d_idae()
    pdae = lotka_volterra_2d_pdae()
    vdae = lotka_volterra_2d_vdae()

    ode_equs = get_function_tuple(ode)
    @test_nowarn ode_equs[:v](ode.t₀, ode.q₀[begin], zero(ode.q₀[begin]))
    @test_nowarn ode_equs[:h](ode.t₀, ode.q₀[begin])

    hode_equs = get_function_tuple(hode)
    @test_nowarn hode_equs[:v](hode.t₀, hode.q₀[begin], hode.p₀[begin], zero(hode.q₀[begin]))
    @test_nowarn hode_equs[:f](hode.t₀, hode.q₀[begin], hode.p₀[begin], zero(hode.p₀[begin]))
    @test_nowarn hode_equs[:h](hode.t₀, hode.q₀[begin], hode.p₀[begin])

    iode_equs = get_function_tuple(iode)
    @test_nowarn iode_equs[:ϑ](iode.t₀, iode.q₀[begin], iode.p₀[begin], zero(iode.q₀[begin]))
    @test_nowarn iode_equs[:f](iode.t₀, iode.q₀[begin], iode.p₀[begin], zero(iode.p₀[begin]))
    @test_nowarn iode_equs[:v̄](iode.t₀, iode.q₀[begin], zero(iode.q₀[begin]))
    @test_nowarn iode_equs[:f̄](iode.t₀, iode.q₀[begin], iode.p₀[begin], zero(iode.p₀[begin]))
    @test_nowarn iode_equs[:h](iode.t₀, iode.q₀[begin])

    pode_equs = get_function_tuple(pode)
    @test_nowarn pode_equs[:v](pode.t₀, pode.q₀[begin], pode.p₀[begin], zero(pode.q₀[begin]))
    @test_nowarn pode_equs[:f](pode.t₀, pode.q₀[begin], pode.p₀[begin], zero(pode.p₀[begin]))
    @test_nowarn pode_equs[:h](pode.t₀, pode.q₀[begin], pode.p₀[begin])

    vode_equs = get_function_tuple(vode)
    @test_nowarn vode_equs[:ϑ](vode.t₀, vode.q₀[begin], vode.p₀[begin], zero(vode.q₀[begin]))
    @test_nowarn vode_equs[:f](vode.t₀, vode.q₀[begin], vode.p₀[begin], zero(vode.p₀[begin]))
    @test_nowarn vode_equs[:v̄](vode.t₀, vode.q₀[begin], zero(vode.q₀[begin]))
    @test_nowarn vode_equs[:f̄](vode.t₀, vode.q₀[begin], vode.p₀[begin], zero(vode.p₀[begin]))
    @test_nowarn vode_equs[:h](vode.t₀, vode.q₀[begin])

    dae_equs = get_function_tuple(dae)
    @test_nowarn dae_equs[:v](dae.t₀, dae.q₀[begin], zero(dae.q₀[begin]))
    @test_nowarn dae_equs[:u](dae.t₀, dae.q₀[begin], dae.λ₀[begin], zero(dae.q₀[begin]))
    @test_nowarn dae_equs[:ϕ](dae.t₀, dae.q₀[begin], zero(dae.λ₀[begin]))
    @test_nowarn dae_equs[:v̄](dae.t₀, dae.q₀[begin], zero(dae.q₀[begin]))
    @test_nowarn dae_equs[:h](dae.t₀, dae.q₀[begin])

    hdae_equs = get_function_tuple(hdae)
    @test_nowarn hdae_equs[:v](hdae.t₀, hdae.q₀[begin], hdae.p₀[begin], zero(hdae.q₀[begin]))
    @test_nowarn hdae_equs[:f](hdae.t₀, hdae.q₀[begin], hdae.p₀[begin], zero(hdae.p₀[begin]))
    @test_nowarn hdae_equs[:u](hdae.t₀, hdae.q₀[begin], hdae.p₀[begin], hdae.λ₀[begin], zero(hdae.q₀[begin]))
    @test_nowarn hdae_equs[:g](hdae.t₀, hdae.q₀[begin], hdae.p₀[begin], hdae.λ₀[begin], zero(hdae.p₀[begin]))
    @test_nowarn hdae_equs[:ϕ](hdae.t₀, hdae.q₀[begin], hdae.p₀[begin], zero(hdae.λ₀[begin]))
    @test_nowarn hdae_equs[:v̄](hdae.t₀, hdae.q₀[begin], hdae.p₀[begin], zero(hdae.q₀[begin]))
    @test_nowarn hdae_equs[:f̄](hdae.t₀, hdae.q₀[begin], hdae.p₀[begin], zero(hdae.p₀[begin]))
    @test_nowarn hdae_equs[:h](hdae.t₀, hdae.q₀[begin], hdae.p₀[begin])

    idae_equs = get_function_tuple(idae)
    @test_nowarn idae_equs[:ϑ](idae.t₀, idae.q₀[begin], idae.p₀[begin], zero(idae.q₀[begin]))
    @test_nowarn idae_equs[:f](idae.t₀, idae.q₀[begin], idae.p₀[begin], zero(idae.p₀[begin]))
    @test_nowarn idae_equs[:u](idae.t₀, idae.q₀[begin], idae.p₀[begin], idae.λ₀[begin], zero(idae.q₀[begin]))
    @test_nowarn idae_equs[:g](idae.t₀, idae.q₀[begin], idae.p₀[begin], idae.λ₀[begin], zero(idae.p₀[begin]))
    @test_nowarn idae_equs[:ϕ](idae.t₀, idae.q₀[begin], idae.p₀[begin], zero(idae.λ₀[begin]))
    @test_nowarn idae_equs[:v̄](idae.t₀, idae.q₀[begin], zero(idae.q₀[begin]))
    @test_nowarn idae_equs[:f̄](idae.t₀, idae.q₀[begin], idae.p₀[begin], zero(idae.p₀[begin]))
    @test_nowarn idae_equs[:h](idae.t₀, idae.q₀[begin])

    pdae_equs = get_function_tuple(pdae)
    @test_nowarn pdae_equs[:v](pdae.t₀, pdae.q₀[begin], pdae.p₀[begin], zero(pdae.q₀[begin]))
    @test_nowarn pdae_equs[:f](pdae.t₀, pdae.q₀[begin], pdae.p₀[begin], zero(pdae.p₀[begin]))
    @test_nowarn pdae_equs[:u](pdae.t₀, pdae.q₀[begin], pdae.p₀[begin], pdae.λ₀[begin], zero(pdae.q₀[begin]))
    @test_nowarn pdae_equs[:g](pdae.t₀, pdae.q₀[begin], pdae.p₀[begin], pdae.λ₀[begin], zero(pdae.p₀[begin]))
    @test_nowarn pdae_equs[:ϕ](pdae.t₀, pdae.q₀[begin], pdae.p₀[begin], zero(pdae.λ₀[begin]))
    @test_nowarn pdae_equs[:v̄](pdae.t₀, pdae.q₀[begin], pdae.p₀[begin], zero(pdae.q₀[begin]))
    @test_nowarn pdae_equs[:f̄](pdae.t₀, pdae.q₀[begin], pdae.p₀[begin], zero(pdae.p₀[begin]))
    @test_nowarn pdae_equs[:h](pdae.t₀, pdae.q₀[begin], pdae.p₀[begin])

    vdae_equs = get_function_tuple(vdae)
    @test_nowarn vdae_equs[:ϑ](vdae.t₀, vdae.q₀[begin], vdae.λ₀[begin], zero(vdae.q₀[begin]))
    @test_nowarn vdae_equs[:f](vdae.t₀, vdae.q₀[begin], vdae.λ₀[begin], zero(vdae.p₀[begin]))
    @test_nowarn vdae_equs[:g](vdae.t₀, vdae.q₀[begin], vdae.λ₀[begin], zero(vdae.p₀[begin]))
    @test_nowarn vdae_equs[:g̅](vdae.t₀, vdae.q₀[begin], vdae.λ₀[begin], zero(vdae.p₀[begin]))
    @test_nowarn vdae_equs[:ϕ](vdae.t₀, vdae.q₀[begin], vdae.p₀[begin], zero(vdae.λ₀[begin]))
    @test_nowarn vdae_equs[:ψ](vdae.t₀, vdae.q₀[begin], vdae.p₀[begin], zero(vdae.q₀[begin]), zero(vdae.p₀[begin]), zero(vdae.λ₀[begin]))
    @test_nowarn vdae_equs[:v̄](vdae.t₀, vdae.q₀[begin], zero(vdae.q₀[begin]))
    @test_nowarn vdae_equs[:f̄](vdae.t₀, vdae.q₀[begin], vdae.λ₀[begin], zero(vdae.p₀[begin]))
    @test_nowarn vdae_equs[:h](vdae.t₀, vdae.q₀[begin])

end


@testset "$(rpad("Kubo Oscillator",80))" begin

    using GeometricProblems.KuboOscillator

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
