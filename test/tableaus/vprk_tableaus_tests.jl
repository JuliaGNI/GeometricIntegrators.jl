
@testset "$(rpad("Variational Partitioned Runge-Kutta Tableaus",80))" begin

    @test typeof(TableauVPGLRK(1)) <: TableauVPRK
    @test typeof(TableauVPGLRK(2)) <: TableauVPRK
    @test typeof(TableauVPGLRK(3)) <: TableauVPRK
    @test typeof(TableauVPGLRK(4)) <: TableauVPRK
    
    @test typeof(TableauVPSRK3()) <: TableauVPRK

    @test typeof(TableauVPLobIIIA(2)) <: TableauVPRK
    @test typeof(TableauVPLobIIIA(3)) <: TableauVPRK
    @test typeof(TableauVPLobIIIA(4)) <: TableauVPRK
    @test typeof(TableauVPLobIIIA(5)) <: TableauVPRK

    @test typeof(TableauVPLobIIIB(2)) <: TableauVPRK
    @test typeof(TableauVPLobIIIB(3)) <: TableauVPRK
    @test typeof(TableauVPLobIIIB(4)) <: TableauVPRK
    @test typeof(TableauVPLobIIIB(5)) <: TableauVPRK

    @test typeof(TableauVPLobIIIC(2)) <: TableauVPRK
    @test typeof(TableauVPLobIIIC(3)) <: TableauVPRK
    @test typeof(TableauVPLobIIIC(4)) <: TableauVPRK
    @test typeof(TableauVPLobIIIC(5)) <: TableauVPRK

    @test typeof(TableauVPLobIIID(2)) <: TableauVPRK
    @test typeof(TableauVPLobIIID(3)) <: TableauVPRK
    @test typeof(TableauVPLobIIID(4)) <: TableauVPRK
    @test typeof(TableauVPLobIIID(5)) <: TableauVPRK

    @test typeof(TableauVPLobIIIE(2)) <: TableauVPRK
    @test typeof(TableauVPLobIIIE(3)) <: TableauVPRK
    @test typeof(TableauVPLobIIIE(4)) <: TableauVPRK
    @test typeof(TableauVPLobIIIE(5)) <: TableauVPRK

    @test typeof(TableauVPLobIIIF(2)) <: TableauVPRK
    @test typeof(TableauVPLobIIIF(3)) <: TableauVPRK
    @test typeof(TableauVPLobIIIF(4)) <: TableauVPRK

    @test typeof(TableauVPLobIIIG(2)) <: TableauVPRK
    @test typeof(TableauVPLobIIIG(3)) <: TableauVPRK
    @test typeof(TableauVPLobIIIG(4)) <: TableauVPRK

    @test typeof(TableauVPLobIIIAIIIA(2)) <: TableauVPRK
    @test typeof(TableauVPLobIIIAIIIA(3)) <: TableauVPRK
    @test typeof(TableauVPLobIIIAIIIA(4)) <: TableauVPRK
    @test typeof(TableauVPLobIIIAIIIA(5)) <: TableauVPRK

    @test typeof(TableauVPLobIIIBIIIB(2)) <: TableauVPRK
    @test typeof(TableauVPLobIIIBIIIB(3)) <: TableauVPRK
    @test typeof(TableauVPLobIIIBIIIB(4)) <: TableauVPRK
    @test typeof(TableauVPLobIIIBIIIB(5)) <: TableauVPRK

    @test typeof(TableauVPRadIIAIIA2()) <: TableauVPRK
    @test typeof(TableauVPRadIIAIIA3()) <: TableauVPRK

end
