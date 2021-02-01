
@testset "$(rpad("Runge-Kutta Tableaus",80))" begin

    # instatiate all partitioned tableaus
    @test typeof(TableauSymplecticEulerA())     <: TableauEPRK
    @test typeof(TableauSymplecticEulerB())     <: TableauEPRK
    @test typeof(TableauLobattoIIIAIIIB2())     <: TableauEPRK
    @test typeof(TableauLobattoIIIBIIIA2())     <: TableauEPRK

    # test instatiation of partioned tableau by composition of two RK tableaus
    @test typeof(TableauEPRK(:PERK4, 4, TableauRK4(), TableauRK4())) <: TableauEPRK
    @test TableauEPRK(:PERK4, 4, TableauRK4(), TableauRK4()) == TableauEPRK(:PERK4, 4, TableauRK4())


    # TODO Add tests for TableauIPRK, TableauSARK, TableauSPARK and TableauGLM.

end
