
@testset "$(rpad("Tableau coefficients",80))" begin

    using GeometricIntegrators.CommonFunctions
    import GeometricIntegrators.Quadratures: LobattoLegendreQuadrature, weights
    import GeometricIntegrators.Tableaus: get_lobatto_nodes, get_lobatto_weights


    # test Gauss-Legendre Runge-Kutta coefficients

    glrk2_tab1 = CoefficientsGLRK(2)
    glrk2_tab2 = CoefficientsGLRK(2, high_precision=false)
    @test glrk2_tab1.a ≈ glrk2_tab2.a atol=4eps()
    @test glrk2_tab1.b == glrk2_tab2.b
    @test glrk2_tab1.c ≈ glrk2_tab2.c atol=4eps()

    glrk3_tab1 = CoefficientsGLRK(3)
    glrk3_tab2 = CoefficientsGLRK(3, high_precision=false)
    @test glrk3_tab1.a ≈ glrk3_tab2.a atol=4eps()
    @test glrk3_tab1.b == glrk3_tab2.b
    @test glrk3_tab1.c ≈ glrk3_tab2.c atol=4eps()

    @test CoefficientsGLRK(1) == getCoefficientsGLRK1()
    @test CoefficientsGLRK(2) == getCoefficientsGLRK2()
    @test CoefficientsGLRK(3) == getCoefficientsGLRK3()
    @test CoefficientsGLRK(4) == getCoefficientsGLRK4()
    @test CoefficientsGLRK(5) == getCoefficientsGLRK5()


    # test Gauss-Lobatto Runge-Kutta coefficients

    function getCoefficientsLobattoIIIA2(T=Float64)
        a = BigFloat[
                [0     0   ]
                [1/2   1/2 ]
            ]

        CoefficientsRK(T, :LobattoIIIA2, 2, a, get_lobatto_weights(2), get_lobatto_nodes(2))
    end

    function getCoefficientsLobattoIIIA3(T=Float64)
        a = BigFloat[
                [0      0     0    ]
                [5/24   1/3  -1/24 ]
                [1/6    2/3   1/6  ]
            ]

        CoefficientsRK(T, :LobattoIIIA3, 4, a, get_lobatto_weights(3), get_lobatto_nodes(3))
    end

    function getCoefficientsLobattoIIIA4(T=Float64)
        a = BigFloat[
                [0            0               0               0          ]
                [(11+√5)/120  (25-   √5)/120  (25-13*√5)/120  (-1+√5)/120]
                [(11-√5)/120  (25+13*√5)/120  (25+   √5)/120  (-1-√5)/120]
                [      1/12            5/12            5/12         1/12 ]
            ]

        CoefficientsRK(T, :LobattoIIIA4, 6, a, get_lobatto_weights(4), get_lobatto_nodes(4))
    end


    function getCoefficientsLobattoIIIA5(T=Float64)
        a = BigFloat[
                [0                 0                   0                  0                   0               ]
                [(119+3*√21)/1960  (343-  9*√21)/2520  (392-96*√21)/2205  (343- 69*√21)/2520  (-21+3*√21)/1960]
                [         13/320   (392+105*√21)/2880             8/45    (392-105*√21)/2880            3/320 ]
                [(119-3*√21)/1960  (343+ 69*√21)/2520  (392+96*√21)/2205  (343+  9*√21)/2520  (-21-3*√21)/1960]
                [          1/20               49/180             16/45               49/180             1/20  ]
            ]

        CoefficientsRK(T, :LobattoIIIA5, 8, a, get_lobatto_weights(5), get_lobatto_nodes(5))
    end


    function getCoefficientsLobattoIIIB2(T=Float64)
        a = BigFloat[
                [1/2  0]
                [1/2  0]
            ]

        CoefficientsRK(T, :LobattoIIIB2, 2, a, get_lobatto_weights(2), get_lobatto_nodes(2))
    end

    function getCoefficientsLobattoIIIB3(T=Float64)
        a = BigFloat[
                [1/6  -1/6   0   ]
                [1/6   1/3   0   ]
                [1/6   5/6   0   ]
            ]

        CoefficientsRK(T, :LobattoIIIB3, 4, a, get_lobatto_weights(3), get_lobatto_nodes(3))
    end

    function getCoefficientsLobattoIIIB4(T=Float64)
        a = BigFloat[
                [ 1/12  (-1-   √5)/24   (-1+   √5)/24    0 ]
                [ 1/12  (25+   √5)/120  (25-13*√5)/120   0 ]
                [ 1/12  (25+13*√5)/120  (25-   √5)/120   0 ]
                [ 1/12  (11-   √5)/24   (11+   √5)/24    0 ]
            ]

        CoefficientsRK(T, :LobattoIIIB4, 6, a, get_lobatto_weights(4), get_lobatto_nodes(4))
    end


    function getCoefficientsLobattoIIIB5(T=Float64)
        a = BigFloat[
                [ 1/20  ( -7-   √21)/120             1/15   ( -7+   √21)/120    0 ]
                [ 1/20  (343+ 9*√21)/2520  (56-15*√21)/315  (343-69*√21)/2520   0 ]
                [ 1/20  ( 49+12*√21)/360             8/45   ( 49-12*√21)/360    0 ]
                [ 1/20  (343+69*√21)/2520  (56+15*√21)/315  (343- 9*√21)/2520   0 ]
                [ 1/20  (119- 3*√21)/360            13/45   (119+ 3*√21)/360    0 ]
            ]

        CoefficientsRK(T, :LobattoIIIB5, 8, a, get_lobatto_weights(5), get_lobatto_nodes(5))
    end


    function getCoefficientsLobattoIIIC2(T=Float64)
        a = BigFloat[
                [1/2  -1/2 ]
                [1/2   1/2 ]
            ]

        CoefficientsRK(T, :LobattoIIIC2, 2, a, get_lobatto_weights(2), get_lobatto_nodes(2))
    end

    function getCoefficientsLobattoIIIC3(T=Float64)
        a = BigFloat[
                [1/6  -1/3    1/6  ]
                [1/6   5/12  -1/12 ]
                [1/6   2/3    1/6  ]
            ]

        CoefficientsRK(T, :LobattoIIIC3, 4, a, get_lobatto_weights(3), get_lobatto_nodes(3))
    end

    function getCoefficientsLobattoIIIC4(T=Float64)
        a = BigFloat[
                [ 1/12        -√5/12         √5/12   -1/12 ]
                [ 1/12          1/4   (10-7*√5)/60   √5/60 ]
                [ 1/12  (10+7*√5)/60          1/4   -√5/60 ]
                [ 1/12          5/12          5/12    1/12 ]
            ]

        CoefficientsRK(T, :LobattoIIIC4, 6, a, get_lobatto_weights(4), get_lobatto_nodes(4))
    end


    function getCoefficientsLobattoIIIC5(T=Float64)
        a = BigFloat[
                [ 1/20            -21/180             2/15             -21/180     1/20  ]
                [ 1/20             29/180   (47-15*√21)/315  (203- 30*√21)/1260   -3/140 ]
                [ 1/20  (329+105*√21)/2880           73/360  (329-105*√21)/2880    3/160 ]
                [ 1/20  (203+ 30*√21)/1260  (47+15*√21)/315             29/180    -3/140 ]
                [ 1/20             49/180            16/45              49/180     1/20  ]
            ]

        CoefficientsRK(T, :LobattoIIIC5, 8, a, get_lobatto_weights(5), get_lobatto_nodes(5))
    end


    function getCoefficientsLobattoIIIC̄2(T=Float64)
        a = BigFloat[
                [0  0]
                [1  0]
            ]

        CoefficientsRK(T, :LobattoIIIC̄2, 2, a, get_lobatto_weights(2), get_lobatto_nodes(2))
    end

    function getCoefficientsLobattoIIIC̄3(T=Float64)
        a = BigFloat[
                [0    0    0 ]
                [1/4  1/4  0 ]
                [0    1    0 ]
            ]

        CoefficientsRK(T, :LobattoIIIC̄3, 4, a, get_lobatto_weights(3), get_lobatto_nodes(3))
    end

    function getCoefficientsLobattoIIIC̄4(T=Float64)
        a = BigFloat[
                [      0             0             0     0 ]
                [ (5+√5)/60          1/6   (15-7*√5)/60  0 ]
                [ (5-√5)/60  (15+7*√5)/60          1/6   0 ]
                [      1/6      (5-√5)/12     (5+√5)/12  0 ]
            ]

        CoefficientsRK(T, :LobattoIIIC̄4, 6, a, get_lobatto_weights(4), get_lobatto_nodes(4))
    end


    function getCoefficientsLobattoIIIC̄5(T=Float64)
        a = BigFloat[
                [ 0     0                0              0                0 ]
                [ 1/14            1/9    (13-3*√21)/63  (14- 3*√21)/126  0 ]
                [ 1/32  (91+21*√21)/576          11/72  (91-21*√21)/576  0 ]
                [ 1/14  (14+ 3*√21)/126  (13+3*√21)/63            1/9    0 ]
                [ 0               7/18            2/9             7/18   0 ]
            ]

        CoefficientsRK(T, :LobattoIIIC̄5, 8, a, get_lobatto_weights(5), get_lobatto_nodes(5))
    end


    function getCoefficientsLobattoIIID2(T=Float64)
        lobC = getCoefficientsLobattoIIIC2(BigFloat)
        lobC̄ = getCoefficientsLobattoIIIC̄2(BigFloat)
        CoefficientsRK(T, :LobattoIIID2, lobC.o, (lobC̄.a + lobC.a)/2, lobC.b, lobC.c)
    end

    function getCoefficientsLobattoIIID3(T=Float64)
        lobC = getCoefficientsLobattoIIIC3(BigFloat)
        lobC̄ = getCoefficientsLobattoIIIC̄3(BigFloat)
        CoefficientsRK(T, :LobattoIIID3, lobC.o, (lobC̄.a + lobC.a)/2, lobC.b, lobC.c)
    end

    function getCoefficientsLobattoIIID4(T=Float64)
        lobC = getCoefficientsLobattoIIIC4(BigFloat)
        lobC̄ = getCoefficientsLobattoIIIC̄4(BigFloat)
        CoefficientsRK(T, :LobattoIIID4, lobC.o, (lobC̄.a + lobC.a)/2, lobC.b, lobC.c)
    end

    function getCoefficientsLobattoIIID5(T=Float64)
        lobC = getCoefficientsLobattoIIIC5(BigFloat)
        lobC̄ = getCoefficientsLobattoIIIC̄5(BigFloat)
        CoefficientsRK(T, :LobattoIIID5, lobC.o, (lobC̄.a + lobC.a)/2, lobC.b, lobC.c)
    end


    function getCoefficientsLobattoIIIE2(T=Float64)
        lobA = getCoefficientsLobattoIIIA2(BigFloat)
        lobB = getCoefficientsLobattoIIIB2(BigFloat)
        CoefficientsRK(T, :LobattoIIIE2, lobA.o, (lobA.a + lobB.a)/2, lobA.b, lobA.c)
    end

    function getCoefficientsLobattoIIIE3(T=Float64)
        lobA = getCoefficientsLobattoIIIA3(BigFloat)
        lobB = getCoefficientsLobattoIIIB3(BigFloat)
        CoefficientsRK(T, :LobattoIIIE3, lobA.o, (lobA.a + lobB.a)/2, lobA.b, lobA.c)
    end

    function getCoefficientsLobattoIIIE4(T=Float64)
        lobA = getCoefficientsLobattoIIIA4(BigFloat)
        lobB = getCoefficientsLobattoIIIB4(BigFloat)
        CoefficientsRK(T, :LobattoIIIE4, lobA.o, (lobA.a + lobB.a)/2, lobA.b, lobA.c)
    end

    function getCoefficientsLobattoIIIE5(T=Float64)
        lobA = getCoefficientsLobattoIIIA5(BigFloat)
        lobB = getCoefficientsLobattoIIIB5(BigFloat)
        CoefficientsRK(T, :LobattoIIIE5, lobA.o, (lobA.a + lobB.a)/2, lobA.b, lobA.c)
    end


    function getCoefficientsLobattoIIIF2(T=Float64)
        a = BigFloat[
                [1/12  -1/12 ]
                [7/12   5/12 ]
            ]

        CoefficientsRK(T, :LobattoIIIF2, 4, a, get_lobatto_weights(2), get_lobatto_nodes(2))
    end

    function getCoefficientsLobattoIIIF3(T=Float64)
        a = BigFloat[
                [1/30   -1/15    1/30 ]
                [5/24    1/3    -1/24 ]
                [2/15   11/15    2/15 ]
            ]

        CoefficientsRK(T, :LobattoIIIF3, 6, a, get_lobatto_weights(3), get_lobatto_nodes(3))
    end

    function getCoefficientsLobattoIIIF4(T=Float64)
        a = BigFloat[
                [  1/56                 -√5/56           √5/56   -1/56         ]
                [ 37/420+√5/120  5/24-   √5/210  5/24-47*√5/420  -1/210+√5/120 ]
                [ 37/420-√5/120  5/24+47*√5/420  5/24+   √5/210  -1/210-√5/120 ]
                [ 17/168         5/12-   √5/56   5/12+   √5/56   11/168        ]
            ]

        CoefficientsRK(T, :LobattoIIIF4, 8, a, get_lobatto_weights(4), get_lobatto_nodes(4))
    end


    function getCoefficientsLobattoIIIG2(T=Float64)
        symplecticize(getCoefficientsLobattoIIIF2(BigFloat); name=:LobattoIIIG2, T=T)
    end

    function getCoefficientsLobattoIIIG3(T=Float64)
        symplecticize(getCoefficientsLobattoIIIF3(BigFloat); name=:LobattoIIIG3, T=T)
    end

    function getCoefficientsLobattoIIIG4(T=Float64)
        symplecticize(getCoefficientsLobattoIIIF4(BigFloat); name=:LobattoIIIG4, T=T)
    end


    @test get_lobatto_nodes(2) ≈ nodes(LobattoLegendreQuadrature(2))
    @test get_lobatto_nodes(3) ≈ nodes(LobattoLegendreQuadrature(3))
    @test get_lobatto_nodes(4) ≈ nodes(LobattoLegendreQuadrature(4))
    @test get_lobatto_nodes(5) ≈ nodes(LobattoLegendreQuadrature(5))
    @test get_lobatto_nodes(6) ≈ nodes(LobattoLegendreQuadrature(6))

    @test get_lobatto_weights(2) ≈ weights(LobattoLegendreQuadrature(2))
    @test get_lobatto_weights(3) ≈ weights(LobattoLegendreQuadrature(3))
    @test get_lobatto_weights(4) ≈ weights(LobattoLegendreQuadrature(4))
    @test get_lobatto_weights(5) ≈ weights(LobattoLegendreQuadrature(5))
    @test get_lobatto_weights(6) ≈ weights(LobattoLegendreQuadrature(6))

    @test CoefficientsLobattoIIIA(2) ≈ getCoefficientsLobattoIIIA2()
    @test CoefficientsLobattoIIIA(3) ≈ getCoefficientsLobattoIIIA3()
    @test CoefficientsLobattoIIIA(4) ≈ getCoefficientsLobattoIIIA4()
    @test CoefficientsLobattoIIIA(5) ≈ getCoefficientsLobattoIIIA5()

    @test CoefficientsLobattoIIIB(2) ≈ getCoefficientsLobattoIIIB2()
    @test CoefficientsLobattoIIIB(3) ≈ getCoefficientsLobattoIIIB3()
    @test CoefficientsLobattoIIIB(4) ≈ getCoefficientsLobattoIIIB4()
    @test CoefficientsLobattoIIIB(5) ≈ getCoefficientsLobattoIIIB5()

    @test CoefficientsLobattoIIIC(2) ≈ getCoefficientsLobattoIIIC2()
    @test CoefficientsLobattoIIIC(3) ≈ getCoefficientsLobattoIIIC3()
    @test CoefficientsLobattoIIIC(4) ≈ getCoefficientsLobattoIIIC4()
    @test CoefficientsLobattoIIIC(5) ≈ getCoefficientsLobattoIIIC5()

    @test CoefficientsLobattoIIIC̄(2) ≈ getCoefficientsLobattoIIIC̄2()
    @test CoefficientsLobattoIIIC̄(3) ≈ getCoefficientsLobattoIIIC̄3()
    @test CoefficientsLobattoIIIC̄(4) ≈ getCoefficientsLobattoIIIC̄4()
    @test CoefficientsLobattoIIIC̄(5) ≈ getCoefficientsLobattoIIIC̄5()

    @test CoefficientsLobattoIIID(2) ≈ getCoefficientsLobattoIIID2()
    @test CoefficientsLobattoIIID(3) ≈ getCoefficientsLobattoIIID3()
    @test CoefficientsLobattoIIID(4) ≈ getCoefficientsLobattoIIID4()
    @test CoefficientsLobattoIIID(5) ≈ getCoefficientsLobattoIIID5()

    @test CoefficientsLobattoIIIE(2) ≈ getCoefficientsLobattoIIIE2()
    @test CoefficientsLobattoIIIE(3) ≈ getCoefficientsLobattoIIIE3()
    @test CoefficientsLobattoIIIE(4) ≈ getCoefficientsLobattoIIIE4()
    @test CoefficientsLobattoIIIE(5) ≈ getCoefficientsLobattoIIIE5()

    @test CoefficientsLobattoIIIF(2) ≈ getCoefficientsLobattoIIIF2()
    @test CoefficientsLobattoIIIF(3) ≈ getCoefficientsLobattoIIIF3()
    @test CoefficientsLobattoIIIF(4) ≈ getCoefficientsLobattoIIIF4()

    @test CoefficientsLobattoIIIG(2) ≈ getCoefficientsLobattoIIIG2()
    @test CoefficientsLobattoIIIG(3) ≈ getCoefficientsLobattoIIIG3()
    @test CoefficientsLobattoIIIG(4) ≈ getCoefficientsLobattoIIIG4()

    
    # test PGLRK coefficients

    @test typeof(CoefficientsPGLRK(2)) <: CoefficientsPGLRK

end
