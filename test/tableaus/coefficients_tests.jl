
@testset "$(rpad("Tableau coefficients",80))" begin

    # test computation of Gauss-Legendre Runge-Kutta coefficients
    glrk2_tab1 = getCoefficientsGLRK(2)
    glrk2_tab2 = getCoefficientsGLRK(2, high_precision=false)
    @test glrk2_tab1.a ≈ glrk2_tab2.a atol=4eps()
    @test glrk2_tab1.b == glrk2_tab2.b
    @test glrk2_tab1.c ≈ glrk2_tab2.c atol=4eps()

    glrk3_tab1 = getCoefficientsGLRK(3)
    glrk3_tab2 = getCoefficientsGLRK(3, high_precision=false)
    @test glrk3_tab1.a ≈ glrk3_tab2.a atol=4eps()
    @test glrk3_tab1.b == glrk3_tab2.b
    @test glrk3_tab1.c ≈ glrk3_tab2.c atol=4eps()

    # test PGLRK coefficients
    pglrk2 = getCoefficientsPGLRK(2)
    @test typeof(pglrk2) <: CoefficientsPGLRK

end
