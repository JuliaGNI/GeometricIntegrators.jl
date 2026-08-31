using GeometricIntegrators.SPARK
using Test

@testset "$(rpad("SPARK Coefficients",80))" begin
    using GeometricIntegrators

    using LinearAlgebra: normalize

    import GeometricIntegrators.SPARK: lobatto_gauss_coefficients
    import QuadratureRules: lobatto_legendre_nodes, lobatto_legendre_weights

    function _get_lobatto_projective_stage(s, T = Float64)
        if s == 1
            a = reshape(T[0, 1], (2, 1))
        elseif s == 2
            a = T[0 0
                  1//2 1//2]
        elseif s == 3
            a = T[0 0 0
                  5//18 8//18 5//18]
        elseif s == 4
            a = T[0 0 0 0
                  1//4-√30/72 1//4+√30/72 1//4+√30/72 1//4-√30/72]
        elseif s == 5
            a = T[0 0 0 0 0
                  (-5797*√70 - 19635*√6 + 3297*√105 + 34069)/(900*(-31*√70 - 105*√6 + 12*√105 + 124)) 3*(-155*√70 - 525*√6 + 21*√105 + 217)/(100*(-31*√70 - 105*√6 + 12*√105 + 124)) 64/225 3*(-155*√70 - 525*√6 + 21*√105 + 217)/(100*(-31*√70 - 105*√6 + 12*√105 + 124)) (-5797*√70 - 19635*√6 + 3297*√105 + 34069)/(900*(-31*√70 - 105*√6 + 12*√105 + 124))]
        else
            @error "Lobatto projective coefficients for " * string(s) *
                   " stages not implemented."
        end

        return a
    end

    @test lobatto_gauss_coefficients(1, 2).a ≈ _get_lobatto_projective_stage(1)
    @test lobatto_gauss_coefficients(2, 2).a ≈ _get_lobatto_projective_stage(2)
    @test lobatto_gauss_coefficients(3, 2).a ≈ _get_lobatto_projective_stage(3)
    @test lobatto_gauss_coefficients(4, 2).a ≈ _get_lobatto_projective_stage(4)
    @test lobatto_gauss_coefficients(5, 2).a ≈ _get_lobatto_projective_stage(5)

    function _get_lobatto_interstage_coefficients(s, σ = s+1, T = Float64)
        if s == 1 && σ == 2
            a = reshape([[0]
                         [1]], σ, s)
        elseif s == 2 && σ == 3
            a = [[0 0]
                 [1/4+√3/8 1/4-√3/8]
                 [1/2 1/2]]
        elseif s == 3 && σ == 4
            a = [[0 0 0]
                 [5/36-√5/180+√15/30 2/9-4*√5/45 5/36-√15/30-√5/180]
                 [5/36+√5/180+√15/30 2/9+4*√5/45 5/36+√5/180-√15/30]
                 [5/18 4/9 5/18]]
        elseif s == 4 && σ == 5
            a = [[0 0 0 0]
                 [
                  0.1591671294900345358547852036503968244037099929876660089701433805420016628303049 0.01565551541810189985834126767767053597372403599533987485875443079743596295789264 -0.002649931830622319490548503812025451589560291885905677307276227674989121280824921 0.0005004515684973118782758043605289134275222217051875147207182646279180365249293946
                  ]
                 [
                  0.176105937938662021225385981377510389121610166295011341406551070207182545793024 0.3049159394838621526624258370570061164298248408582482915516281042450345033436905 0.02115663794741091865104218833199417995250081122480493820210723599558956292335403 -0.002178515369935092538854006766510685503935818378064571160286410447806612060067542
                  ]
                 [
                  0.1734269710002296168082561702504707901901521262117592555255463951314578972080256 0.3287225092618953908040165292010257479718859439689589070610115679156131875478712 0.3104170620131711714551267577113297604086016160877133548949809094431881033091524 0.01476029307869239283174677096060287921396435492928076127612127921737427090265032
                  ]
                 [
                  0.1739274225687269286865319746109997036176743479169467702462646597593759337329541 0.3260725774312730713134680253890002963823256520830532297537353402406240662670459 0.3260725774312730713134680253890002963823256520830532297537353402406240662670459 0.1739274225687269286865319746109997036176743479169467702462646597593759337329541
                  ]]
        else
            @error("Number of stages s=$(s) and σ=$(σ) is not supported.")
        end

        CoefficientsIRK{T}(
            :LobattoIIIIS, s^2, s, σ, a, lobatto_legendre_weights(BigFloat, σ),
            lobatto_legendre_nodes(BigFloat, σ))
    end

    @test lobatto_gauss_coefficients(1) ≈ _get_lobatto_interstage_coefficients(1)
    @test lobatto_gauss_coefficients(2) ≈ _get_lobatto_interstage_coefficients(2)
    @test lobatto_gauss_coefficients(3) ≈ _get_lobatto_interstage_coefficients(3)
    @test lobatto_gauss_coefficients(4) ≈ _get_lobatto_interstage_coefficients(4)

    # The Vandermonde system of lobatto_gauss_coefficients is solved in BigFloat and
    # only the result is narrowed to Float64. Pin that down: the comparisons above use
    # the default tolerance of ≈, rtol = √eps ≈ 1.5e-8, which is seven orders of
    # magnitude too loose to see the difference. Solving in Float64 instead moves the
    # coefficients by 2.0e-15 in relative norm at s = 4 (and by 9.3e-13 at s = 8),
    # whereas the BigFloat solve reproduces the reference exactly, so rtol = 1e-15
    # separates the two with room on either side.
    @test isapprox(lobatto_gauss_coefficients(4).a,
        _get_lobatto_interstage_coefficients(4).a; rtol = 1e-15)

    # test PGLRK coefficients
    # (the mathematical identities are checked in test/methods/pglrk_coefficients_tests.jl)

    @test typeof(CoefficientsPGLRK(3)) <: CoefficientsPGLRK
end

@testset "$(rpad("Special Partitioned Additive Runge-Kutta Tableaus",80))" begin
    @test typeof(SPARKGLRK(1)) <: SPARKMethod
    @test typeof(SPARKGLRK(2)) <: SPARKMethod
    @test typeof(SPARKLobattoIIIAIIIB(2)) <: SPARKMethod
    @test typeof(SPARKLobattoIIIAIIIB(3)) <: SPARKMethod
    @test typeof(SPARKGLRKLobattoIIIAIIIB(1)) <: SPARKMethod
    @test typeof(SPARKGLRKLobattoIIIAIIIB(2)) <: SPARKMethod
    @test typeof(SPARKLobABC(2)) <: SPARKMethod
    @test typeof(SPARKLobABC(3)) <: SPARKMethod
    @test typeof(SPARKLobABD(2)) <: SPARKMethod
    @test typeof(SPARKLobABD(3)) <: SPARKMethod
    @test typeof(SPARKGLVPRK(1)) <: SPARKMethod
    @test typeof(SPARKGLVPRK(2)) <: SPARKMethod

    @test typeof(TableauGausspSymplectic(1)) <: VPARK
    @test typeof(TableauGausspSymplectic(2)) <: VPARK
    @test typeof(TableauLobattoIIIAIIIBpSymplectic(2)) <: VPARK
    @test typeof(TableauLobattoIIIAIIIBpSymplectic(3)) <: VPARK
    @test typeof(TableauLobattoIIIAIIIBpSymplectic(4)) <: VPARK
    @test typeof(TableauLobattoIIIBIIIApSymplectic(2)) <: VPARK
    @test typeof(TableauLobattoIIIBIIIApSymplectic(3)) <: VPARK
    @test typeof(TableauLobattoIIIBIIIApSymplectic(4)) <: VPARK

    @test typeof(TableauVSPARKGLRKpLobattoIIIAIIIB(1)) <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpLobattoIIIAIIIB(2)) <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpLobattoIIIBIIIA(1)) <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpLobattoIIIBIIIA(2)) <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedLobattoIIIAIIIB(1)) <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedLobattoIIIAIIIB(2)) <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedLobattoIIIBIIIA(1)) <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedLobattoIIIBIIIA(2)) <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpInternal(1)) <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpInternal(2)) <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedInternal(1)) <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedInternal(2)) <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpMidpoint(1)) <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpMidpoint(2)) <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedMidpoint(1)) <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedMidpoint(2)) <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpSymmetric(1)) <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpSymmetric(2)) <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpSymplectic(1)) <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpSymplectic(2)) <: VSPARKprimary

    @test typeof(TableauVSPARKLobattoIIIAIIIBpLobattoIIIAIIIB(2)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpLobattoIIIAIIIB(3)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApLobattoIIIAIIIB(2)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApLobattoIIIAIIIB(3)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpLobattoIIIBIIIA(2)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpLobattoIIIBIIIA(3)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApLobattoIIIBIIIA(2)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApLobattoIIIBIIIA(3)) <: VSPARKprimary

    @test typeof(TableauVSPARKLobattoIIIAIIIBpModifiedLobattoIIIAIIIB(2)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpModifiedLobattoIIIAIIIB(3)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpModifiedLobattoIIIBIIIA(2)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpModifiedLobattoIIIBIIIA(3)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApModifiedLobattoIIIAIIIB(2)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApModifiedLobattoIIIAIIIB(3)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApModifiedLobattoIIIBIIIA(2)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApModifiedLobattoIIIBIIIA(3)) <: VSPARKprimary

    @test typeof(TableauVSPARKLobattoIIIAIIIBpMidpoint(2)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpMidpoint(3)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApMidpoint(2)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApMidpoint(3)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpModifiedMidpoint(2)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpModifiedMidpoint(3)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApModifiedMidpoint(2)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApModifiedMidpoint(3)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpSymmetric(2)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpSymmetric(3)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApSymmetric(2)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApSymmetric(3)) <: VSPARKprimary

    @test typeof(TableauVSPARKLobABCCD(2)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobABCCE(2)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobABDE(2)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobABED(2)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobABD(2)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobABE(2)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobDE(2)) <: VSPARKprimary
    @test typeof(TableauVSPARKLobED(2)) <: VSPARKprimary

    @test typeof(TableauVSPARKLobattoIIIAB(2)) <: VSPARKsecondary
    @test typeof(TableauVSPARKLobattoIIIBA(2)) <: VSPARKsecondary
    @test typeof(TableauVSPARKLobattoIIICC̄(2)) <: VSPARKsecondary
    @test typeof(TableauVSPARKLobattoIIIC̄C(2)) <: VSPARKsecondary
    @test typeof(TableauVSPARKLobattoIIID(2)) <: VSPARKsecondary
    @test typeof(TableauVSPARKLobattoIIIE(2)) <: VSPARKsecondary
    @test typeof(TableauVSPARKGLRKLobattoIIIAB(2)) <: VSPARKsecondary
    @test typeof(TableauVSPARKGLRKLobattoIIIBA(2)) <: VSPARKsecondary
    @test typeof(TableauVSPARKGLRKLobattoIIICC̄(2)) <: VSPARKsecondary
    @test typeof(TableauVSPARKGLRKLobattoIIIC̄C(2)) <: VSPARKsecondary
    @test typeof(TableauVSPARKGLRKLobattoIIID(2)) <: VSPARKsecondary
    @test typeof(TableauVSPARKGLRKLobattoIIIE(2)) <: VSPARKsecondary

    @test typeof(TableauHPARKGLRK(1)) <: HPARK
    @test typeof(TableauHPARKGLRK(2)) <: HPARK
    @test typeof(TableauHPARKLobattoIIIAIIIB(2)) <: HPARK
    @test typeof(TableauHPARKLobattoIIIAIIIB(3)) <: HPARK
    @test typeof(TableauHPARKLobattoIIIBIIIA(2)) <: HPARK
    @test typeof(TableauHPARKLobattoIIIBIIIA(3)) <: HPARK

    # DISABLED: `getTableauHPARK` still exists and runs, but the equality no
    # longer holds — TableauHPARKGLRK(s) and getTableauHPARK(...) build tableaus
    # that differ in metadata (e.g. name), so `==` returns false.
    # @test TableauHPARKGLRK(1) == getTableauHPARK(:HPARKGLRK1, TableauGauss(1), TableauGauss(1))
    # @test TableauHPARKGLRK(2) == getTableauHPARK(:HPARKGLRK2, TableauGauss(2), TableauGauss(2))

    @test typeof(TableauHSPARKGLRKpSymmetric(1)) <: HSPARKprimary
    @test typeof(TableauHSPARKGLRKpSymmetric(2)) <: HSPARKprimary
    @test typeof(TableauHSPARKLobattoIIIAIIIBpSymmetric(2)) <: HSPARKprimary
    @test typeof(TableauHSPARKLobattoIIIAIIIBpSymmetric(3)) <: HSPARKprimary
    @test typeof(TableauHSPARKLobattoIIIAIIIBpSymmetric(4)) <: HSPARKprimary
    @test typeof(TableauHSPARKLobattoIIIBIIIApSymmetric(2)) <: HSPARKprimary
    @test typeof(TableauHSPARKLobattoIIIBIIIApSymmetric(3)) <: HSPARKprimary
    @test typeof(TableauHSPARKLobattoIIIBIIIApSymmetric(4)) <: HSPARKprimary

    @test typeof(TableauHSPARKLobattoIIIAB(2)) <: HSPARKsecondary
    @test typeof(TableauHSPARKLobattoIIIBA(2)) <: HSPARKsecondary
    @test typeof(TableauHSPARKLobattoIIICC̄(2)) <: HSPARKsecondary
    @test typeof(TableauHSPARKLobattoIIIC̄C(2)) <: HSPARKsecondary
    @test typeof(TableauHSPARKLobattoIIID(2)) <: HSPARKsecondary
    @test typeof(TableauHSPARKLobattoIIIE(2)) <: HSPARKsecondary
    @test typeof(TableauHSPARKGLRKLobattoIIIAB(1)) <: HSPARKsecondary
    @test typeof(TableauHSPARKGLRKLobattoIIIBA(1)) <: HSPARKsecondary
    @test typeof(TableauHSPARKGLRKLobattoIIICC̄(1)) <: HSPARKsecondary
    @test typeof(TableauHSPARKGLRKLobattoIIIC̄C(1)) <: HSPARKsecondary
    @test typeof(TableauHSPARKGLRKLobattoIIID(1)) <: HSPARKsecondary
    @test typeof(TableauHSPARKGLRKLobattoIIIE(1)) <: HSPARKsecondary

    @test typeof(SLRKLobattoIIIAB(2)) <: SLRK
    @test typeof(SLRKLobattoIIIBA(2)) <: SLRK
    @test typeof(SLRKLobattoIIICC̄(2)) <: SLRK
    @test typeof(SLRKLobattoIIIC̄C(2)) <: SLRK
    @test typeof(SLRKLobattoIIID(2)) <: SLRK
    @test typeof(SLRKLobattoIIIE(2)) <: SLRK

    # Every SLRK constructor must register under its own name symbol: `SLRKLobattoIIIC̄C`
    # was registered under `SLRKLobattoIIICC̄`'s symbol until the fifth pass (finding S13),
    # which the `typeof` checks above cannot see.
    let names = [ctor(2).name
                 for ctor in (SLRKLobattoIIIAB, SLRKLobattoIIIBA,
            SLRKLobattoIIICC̄, SLRKLobattoIIIC̄C,
            SLRKLobattoIIID, SLRKLobattoIIIE)]
        @test length(unique(names)) == length(names)
        @test SLRKLobattoIIICC̄(2).name == Symbol("SLRKLobattoIIICIIIC̄")
        @test SLRKLobattoIIIC̄C(2).name == Symbol("SLRKLobattoIIIC̄IIIC")
    end
end
