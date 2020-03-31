
@testset "$(rpad("Special Partitioned Additive Runge-Kutta Tableaus",80))" begin

    @test typeof(getTableauGLRKpSymplectic(1)) <: TableauVPARK
    @test typeof(getTableauGLRKpSymplectic(2)) <: TableauVPARK
    @test typeof(getTableauLobIIIAIIIB2pSymplectic()) <: TableauVPARK
    @test typeof(getTableauLobIIIAIIIB3pSymplectic()) <: TableauVPARK
    @test typeof(getTableauLobIIIAIIIB4pSymplectic()) <: TableauVPARK

    @test typeof(getTableauVSPARKGLRKpMidpoint(1))   <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKGLRKpMidpoint(2))   <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKGLRKpSymmetric(1))  <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKGLRKpSymmetric(2))  <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKGLRKpSymplectic(1)) <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKGLRKpSymplectic(1)) <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobIIIAIIIB2pSymmetric()) <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobIIIAIIIB3pSymmetric()) <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobIIIAIIIB4pSymmetric()) <: TableauVSPARKprimary

    @test typeof(getTableauVSPARKLobIIIAB(2)) <: TableauVSPARKsecondary
    @test typeof(getTableauVSPARKLobIIIC(2))  <: TableauVSPARKsecondary
    @test typeof(getTableauVSPARKLobIIID(2))  <: TableauVSPARKsecondary
    @test typeof(getTableauVSPARKLobIIIE(2))  <: TableauVSPARKsecondary
    @test typeof(getTableauVSPARKGLRKLobIIIAB(2)) <: TableauVSPARKsecondary
    @test typeof(getTableauVSPARKGLRKLobIIIC(2))  <: TableauVSPARKsecondary
    @test typeof(getTableauVSPARKGLRKLobIIID(2))  <: TableauVSPARKsecondary
    @test typeof(getTableauVSPARKGLRKLobIIIE(2))  <: TableauVSPARKsecondary

    @test typeof(getTableauHPARKGLRK(1)) <: TableauHPARK
    @test typeof(getTableauHPARKGLRK(2)) <: TableauHPARK
    @test typeof(getTableauHPARKLobIIIAIIIB2()) <: TableauHPARK
    @test typeof(getTableauHPARKLobIIIAIIIB3()) <: TableauHPARK
    @test typeof(getTableauHPARKLobIIIAIIIB4()) <: TableauHPARK

    @test typeof(getTableauHSPARKGLRKpSymmetric(1)) <: TableauHSPARKprimary
    @test typeof(getTableauHSPARKGLRKpSymmetric(2)) <: TableauHSPARKprimary
    @test typeof(getTableauHSPARKLobIIIAIIIB2pSymmetric()) <: TableauHSPARKprimary
    @test typeof(getTableauHSPARKLobIIIAIIIB3pSymmetric()) <: TableauHSPARKprimary
    @test typeof(getTableauHSPARKLobIIIAIIIB4pSymmetric()) <: TableauHSPARKprimary

    @test typeof(getTableauHSPARKLobIIIAB(2)) <: TableauHSPARKsecondary
    @test typeof(getTableauHSPARKLobIIIC(2))  <: TableauHSPARKsecondary
    @test typeof(getTableauHSPARKLobIIID(2))  <: TableauHSPARKsecondary
    @test typeof(getTableauHSPARKLobIIIE(2))  <: TableauHSPARKsecondary
    @test typeof(getTableauHSPARKGLRKLobIIIAB(1)) <: TableauHSPARKsecondary
    @test typeof(getTableauHSPARKGLRKLobIIIC(1))  <: TableauHSPARKsecondary
    @test typeof(getTableauHSPARKGLRKLobIIID(1))  <: TableauHSPARKsecondary
    @test typeof(getTableauHSPARKGLRKLobIIIE(1))  <: TableauHSPARKsecondary

    @test typeof(getTableauSLRKLobIIIAB(2)) <: TableauSLRK
    @test typeof(getTableauSLRKLobIIIC(2))  <: TableauSLRK
    @test typeof(getTableauSLRKLobIIID(2))  <: TableauSLRK
    @test typeof(getTableauSLRKLobIIIE(2))  <: TableauSLRK

end
