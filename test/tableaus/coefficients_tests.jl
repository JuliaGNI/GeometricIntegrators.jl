
@testset "$(rpad("Tableau coefficients",80))" begin

    using LinearAlgebra: normalize
    using GeometricIntegrators.CommonFunctions
    import GeometricIntegrators.Quadratures: LobattoLegendreQuadrature, weights
    import GeometricIntegrators.Tableaus: get_lobatto_nodes, get_lobatto_weights,
                                          get_lobatto_glrk_coefficients, get_lobatto_d_vector


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


    function _get_lobatto_projective_stage(s, T=Float64)
        if s == 1
            a = reshape(T[0,1], (2,1))
        elseif s == 2
            a = T[ 0     0
                1//2  1//2 ]
        elseif s == 3
            a = T[ 0      0      0
                5//18  8//18  5//18 ]
        elseif s == 4
            a = T[ 0            0            0            0
                1//4-√30/72  1//4+√30/72  1//4+√30/72  1//4-√30/72 ]
        elseif s == 5
            a = T[ 0            0            0            0            0
                (-5797*√70 - 19635*√6 + 3297*√105 + 34069)/(900*(-31*√70 - 105*√6 + 12*√105 + 124))  3*(-155*√70 - 525*√6 + 21*√105 + 217)/(100*(-31*√70 - 105*√6 + 12*√105 + 124))  64/225  3*(-155*√70 - 525*√6 + 21*√105 + 217)/(100*(-31*√70 - 105*√6 + 12*√105 + 124))  (-5797*√70 - 19635*√6 + 3297*√105 + 34069)/(900*(-31*√70 - 105*√6 + 12*√105 + 124)) ]
        else
            @error "Lobatto projective coefficients for " * string(s) * " stages not implemented."
        end
    
        return a
    end

    @test get_lobatto_glrk_coefficients(1,2).a ≈ _get_lobatto_projective_stage(1)
    @test get_lobatto_glrk_coefficients(2,2).a ≈ _get_lobatto_projective_stage(2)
    @test get_lobatto_glrk_coefficients(3,2).a ≈ _get_lobatto_projective_stage(3)
    @test get_lobatto_glrk_coefficients(4,2).a ≈ _get_lobatto_projective_stage(4)
    @test get_lobatto_glrk_coefficients(5,2).a ≈ _get_lobatto_projective_stage(5)


    function _get_lobatto_interstage_coefficients(s, σ=s+1, T=Float64)
        if s == 1 && σ == 2
            a = reshape([
                    [0]
                    [1]
                ], σ, s)
        elseif s == 2 && σ == 3
            a = [
                    [0         0       ]
                    [1/4+√3/8  1/4-√3/8]
                    [1/2       1/2     ]
                ]
        elseif s == 3 && σ == 4
            a = [
                    [0                     0              0                 ]
                    [5/36-√5/180+√15/30    2/9-4*√5/45    5/36-√15/30-√5/180]
                    [5/36+√5/180+√15/30    2/9+4*√5/45    5/36+√5/180-√15/30]
                    [5/18                  4/9            5/18              ]
                ]
        elseif s == 4 && σ == 5
            a = [
                    [ 0  0  0  0 ]
                    [ 0.1591671294900345358547852036503968244037099929876660089701433805420016628303049  0.01565551541810189985834126767767053597372403599533987485875443079743596295789264  -0.002649931830622319490548503812025451589560291885905677307276227674989121280824921  0.0005004515684973118782758043605289134275222217051875147207182646279180365249293946 ]
                    [ 0.176105937938662021225385981377510389121610166295011341406551070207182545793024   0.3049159394838621526624258370570061164298248408582482915516281042450345033436905    0.02115663794741091865104218833199417995250081122480493820210723599558956292335403  -0.002178515369935092538854006766510685503935818378064571160286410447806612060067542  ]
                    [ 0.1734269710002296168082561702504707901901521262117592555255463951314578972080256  0.3287225092618953908040165292010257479718859439689589070610115679156131875478712    0.3104170620131711714551267577113297604086016160877133548949809094431881033091524    0.01476029307869239283174677096060287921396435492928076127612127921737427090265032   ]
                    [ 0.1739274225687269286865319746109997036176743479169467702462646597593759337329541  0.3260725774312730713134680253890002963823256520830532297537353402406240662670459    0.3260725774312730713134680253890002963823256520830532297537353402406240662670459    0.1739274225687269286865319746109997036176743479169467702462646597593759337329541    ]
                ]            
        else
            @error("Number of stages s=$(s) and σ=$(σ) is not supported.")
        end
 
        CoefficientsIRK{T}(:LobattoIIIIS, s^2, s, σ, a, get_lobatto_weights(σ), get_lobatto_nodes(σ))
    end

    @test get_lobatto_glrk_coefficients(1) ≈ _get_lobatto_interstage_coefficients(1)
    @test get_lobatto_glrk_coefficients(2) ≈ _get_lobatto_interstage_coefficients(2)
    @test get_lobatto_glrk_coefficients(3) ≈ _get_lobatto_interstage_coefficients(3)
    @test get_lobatto_glrk_coefficients(4) ≈ _get_lobatto_interstage_coefficients(4)


    function _get_lobatto_d_vector(s)
        if s == 2
            d = [+1.0, -1.0]
        elseif s == 3
            d = [+1.0, -2.0, +1.0]
        elseif s == 4
            d = [+1.0, -√5, +√5, -1.0]
        elseif s == 5
            d = [+3.0, -7.0, +8.0, -7.0, +3.0]
        else
            @error("We don't have a d vector for s=$(s) stages.")
        end
        return normalize(d) * sign(d[begin])
    end

    @test get_lobatto_d_vector(2; normalize=true) ≈ _get_lobatto_d_vector(2)
    @test get_lobatto_d_vector(3; normalize=true) ≈ _get_lobatto_d_vector(3)
    @test get_lobatto_d_vector(4; normalize=true) ≈ _get_lobatto_d_vector(4)
    @test get_lobatto_d_vector(5; normalize=true) ≈ _get_lobatto_d_vector(5)


    # test PGLRK coefficients

    @test typeof(CoefficientsPGLRK(2)) <: CoefficientsPGLRK

end
