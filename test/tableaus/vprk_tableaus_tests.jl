
@testset "$(rpad("Variational Partitioned Runge-Kutta Tableaus",80))" begin

    @test typeof(getTableauVPGLRK(1)) <: TableauVPRK
    @test typeof(getTableauVPGLRK(2)) <: TableauVPRK
    @test typeof(getTableauVPSRK3()) <: TableauVPRK

    @test typeof(getTableauVPLobIIIA2()) <: TableauVPRK
    @test typeof(getTableauVPLobIIIA3()) <: TableauVPRK
    @test typeof(getTableauVPLobIIIA4()) <: TableauVPRK
    @test typeof(getTableauVPLobIIIB2()) <: TableauVPRK
    @test typeof(getTableauVPLobIIIB3()) <: TableauVPRK
    @test typeof(getTableauVPLobIIIB4()) <: TableauVPRK
    @test typeof(getTableauVPLobIIIC2()) <: TableauVPRK
    @test typeof(getTableauVPLobIIIC3()) <: TableauVPRK
    @test typeof(getTableauVPLobIIIC4()) <: TableauVPRK
    @test typeof(getTableauVPLobIIID2()) <: TableauVPRK
    @test typeof(getTableauVPLobIIID3()) <: TableauVPRK
    @test typeof(getTableauVPLobIIID4()) <: TableauVPRK
    @test typeof(getTableauVPLobIIIE2()) <: TableauVPRK
    @test typeof(getTableauVPLobIIIE3()) <: TableauVPRK
    @test typeof(getTableauVPLobIIIE4()) <: TableauVPRK
    @test typeof(getTableauVPLobIIIF2()) <: TableauVPRK
    @test typeof(getTableauVPLobIIIF3()) <: TableauVPRK
    @test typeof(getTableauVPLobIIIF4()) <: TableauVPRK
    @test typeof(getTableauVPLobIIIG2()) <: TableauVPRK
    @test typeof(getTableauVPLobIIIG3()) <: TableauVPRK
    @test typeof(getTableauVPLobIIIG4()) <: TableauVPRK

    @test typeof(getTableauVPLobIIIAIIIA2()) <: TableauVPRK
    @test typeof(getTableauVPLobIIIAIIIA3()) <: TableauVPRK
    @test typeof(getTableauVPLobIIIAIIIA4()) <: TableauVPRK

    @test typeof(getTableauVPRadIIAIIA2()) <: TableauVPRK
    @test typeof(getTableauVPRadIIAIIA3()) <: TableauVPRK

end
