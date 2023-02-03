
@testset "$(rpad("Tableau coefficients",80))" begin

    using LinearAlgebra: normalize
    using GeometricBase
    using GeometricIntegrators.Utils
    
    import GeometricIntegrators.Tableaus: get_lobatto_nodes, get_lobatto_weights,
                                          lobatto_gauss_coefficients
    import QuadratureRules: LobattoLegendreQuadrature, nodes, weights



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

    @test lobatto_gauss_coefficients(1,2).a ≈ _get_lobatto_projective_stage(1)
    @test lobatto_gauss_coefficients(2,2).a ≈ _get_lobatto_projective_stage(2)
    @test lobatto_gauss_coefficients(3,2).a ≈ _get_lobatto_projective_stage(3)
    @test lobatto_gauss_coefficients(4,2).a ≈ _get_lobatto_projective_stage(4)
    @test lobatto_gauss_coefficients(5,2).a ≈ _get_lobatto_projective_stage(5)


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

    @test lobatto_gauss_coefficients(1) ≈ _get_lobatto_interstage_coefficients(1)
    @test lobatto_gauss_coefficients(2) ≈ _get_lobatto_interstage_coefficients(2)
    @test lobatto_gauss_coefficients(3) ≈ _get_lobatto_interstage_coefficients(3)
    @test lobatto_gauss_coefficients(4) ≈ _get_lobatto_interstage_coefficients(4)


    # test PGLRK coefficients

    # TODO: reactivate
    # @test typeof(CoefficientsPGLRK(2)) <: CoefficientsPGLRK

end
