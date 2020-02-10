
@testset "$(rpad("Tableau coefficients",80))" begin

    # test Gauss-Legendre Runge-Kutta coefficients

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

    @test getCoefficientsGLRK(1) == getCoefficientsGLRK1()
    @test getCoefficientsGLRK(2) == getCoefficientsGLRK2()
    @test getCoefficientsGLRK(3) == getCoefficientsGLRK3()
    @test getCoefficientsGLRK(4) == getCoefficientsGLRK4()
    @test getCoefficientsGLRK(5) == getCoefficientsGLRK5()


    # test Gauss-Lobatto Runge-Kutta coefficients

    @test getCoefficientsLobIII(2) == getCoefficientsLobIII2()
    @test getCoefficientsLobIII(3) == getCoefficientsLobIII3()
    @test getCoefficientsLobIII(4) == getCoefficientsLobIII4()

    @test getCoefficientsLobIIIA(2) == getCoefficientsLobIIIA2()
    @test getCoefficientsLobIIIA(3) == getCoefficientsLobIIIA3()
    @test getCoefficientsLobIIIA(4) == getCoefficientsLobIIIA4()

    @test getCoefficientsLobIIIB(2) == getCoefficientsLobIIIB2()
    @test getCoefficientsLobIIIB(3) == getCoefficientsLobIIIB3()
    @test getCoefficientsLobIIIB(4) == getCoefficientsLobIIIB4()

    @test getCoefficientsLobIIIC(2) == getCoefficientsLobIIIC2()
    @test getCoefficientsLobIIIC(3) == getCoefficientsLobIIIC3()
    @test getCoefficientsLobIIIC(4) == getCoefficientsLobIIIC4()

    @test getCoefficientsLobIIID(2) == getCoefficientsLobIIID2()
    @test getCoefficientsLobIIID(3) == getCoefficientsLobIIID3()
    @test getCoefficientsLobIIID(4) == getCoefficientsLobIIID4()

    @test getCoefficientsLobIIIE(2) == getCoefficientsLobIIIE2()
    @test getCoefficientsLobIIIE(3) == getCoefficientsLobIIIE3()
    @test getCoefficientsLobIIIE(4) == getCoefficientsLobIIIE4()


    # test PGLRK coefficients

    pglrk2 = getCoefficientsPGLRK(2)
    @test typeof(pglrk2) <: CoefficientsPGLRK

end
