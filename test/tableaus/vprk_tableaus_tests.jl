
@testset "$(rpad("Variational Partitioned Runge-Kutta Tableaus",80))" begin

    @test typeof(TableauVPGLRK(1)) <: TableauVPRK
    @test typeof(TableauVPGLRK(2)) <: TableauVPRK
    @test typeof(TableauVPGLRK(3)) <: TableauVPRK
    @test typeof(TableauVPGLRK(4)) <: TableauVPRK

    @test typeof(TableauVPSRK3()) <: TableauVPRK

    @test typeof(TableauVPRadauIIA(2)) <: TableauVPRK
    @test typeof(TableauVPRadauIIA(3)) <: TableauVPRK
    @test typeof(TableauVPRadauIIA(4)) <: TableauVPRK
    @test typeof(TableauVPRadauIIA(5)) <: TableauVPRK

    @test typeof(TableauVPRadauIIB(2)) <: TableauVPRK
    @test typeof(TableauVPRadauIIB(3)) <: TableauVPRK
    @test typeof(TableauVPRadauIIB(4)) <: TableauVPRK
    @test typeof(TableauVPRadauIIB(5)) <: TableauVPRK

    @test typeof(TableauVPLobattoIIIA(2)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIA(3)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIA(4)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIA(5)) <: TableauVPRK

    @test typeof(TableauVPLobattoIIIB(2)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIB(3)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIB(4)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIB(5)) <: TableauVPRK

    @test typeof(TableauVPLobattoIIIC(2)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIC(3)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIC(4)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIC(5)) <: TableauVPRK

    @test typeof(TableauVPLobattoIIIC̄(2)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIC̄(3)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIC̄(4)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIC̄(5)) <: TableauVPRK

    @test typeof(TableauVPLobattoIIID(2)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIID(3)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIID(4)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIID(5)) <: TableauVPRK

    @test typeof(TableauVPLobattoIIIE(2)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIE(3)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIE(4)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIE(5)) <: TableauVPRK

    @test typeof(TableauVPLobattoIIIF(2)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIF(3)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIF(4)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIF(5)) <: TableauVPRK

    @test typeof(TableauVPLobattoIIIG(2)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIG(3)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIG(4)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIG(5)) <: TableauVPRK

    @test typeof(TableauVPLobattoIIIAIIIB(2)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIAIIIB(3)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIAIIIB(4)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIAIIIB(5)) <: TableauVPRK

    @test typeof(TableauVPLobattoIIIBIIIA(2)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIBIIIA(3)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIBIIIA(4)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIBIIIA(5)) <: TableauVPRK

    @test typeof(TableauVPLobattoIIIAIIIĀ(2)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIAIIIĀ(3)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIAIIIĀ(4)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIAIIIĀ(5)) <: TableauVPRK

    @test typeof(TableauVPLobattoIIIBIIIB̄(2)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIBIIIB̄(3)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIBIIIB̄(4)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIBIIIB̄(5)) <: TableauVPRK

    @test typeof(TableauVPLobattoIIICIIIC̄(2)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIICIIIC̄(3)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIICIIIC̄(4)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIICIIIC̄(5)) <: TableauVPRK

    @test typeof(TableauVPLobattoIIIDIIID̄(2)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIDIIID̄(3)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIDIIID̄(4)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIDIIID̄(5)) <: TableauVPRK

    @test typeof(TableauVPLobattoIIIEIIIĒ(2)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIEIIIĒ(3)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIEIIIĒ(4)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIEIIIĒ(5)) <: TableauVPRK

    @test typeof(TableauVPLobattoIIIFIIIF̄(2)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIFIIIF̄(3)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIFIIIF̄(4)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIFIIIF̄(5)) <: TableauVPRK

    @test typeof(TableauVPLobattoIIIGIIIḠ(2)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIGIIIḠ(3)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIGIIIḠ(4)) <: TableauVPRK
    @test typeof(TableauVPLobattoIIIGIIIḠ(5)) <: TableauVPRK


    @test TableauVPSRK3()         == TableauVPSRK3(Float64)
    @test TableauVPGLRK(1)        == TableauVPGLRK(Float64, 1)
    @test TableauVPRadauIIA(2)    == TableauVPRadauIIA(Float64, 2)
    @test TableauVPRadauIIB(2)    == TableauVPRadauIIB(Float64, 2)
    @test TableauVPLobattoIIIA(2) == TableauVPLobattoIIIA(Float64, 2)
    @test TableauVPLobattoIIIB(2) == TableauVPLobattoIIIB(Float64, 2)
    @test TableauVPLobattoIIIC(2) == TableauVPLobattoIIIC(Float64, 2)
    @test TableauVPLobattoIIIC̄(2) == TableauVPLobattoIIIC̄(Float64, 2)
    @test TableauVPLobattoIIID(2) == TableauVPLobattoIIID(Float64, 2)
    @test TableauVPLobattoIIIE(2) == TableauVPLobattoIIIE(Float64, 2)
    @test TableauVPLobattoIIIF(2) == TableauVPLobattoIIIF(Float64, 2)
    @test TableauVPLobattoIIIG(2) == TableauVPLobattoIIIG(Float64, 2)

end
