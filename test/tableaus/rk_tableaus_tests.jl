
@testset "$(rpad("Runge-Kutta Tableaus",80))" begin

    # instatiate all partitioned tableaus
    @test typeof(TableauSymplecticEulerA())     <: TableauPRK
    @test typeof(TableauSymplecticEulerB())     <: TableauPRK
    @test typeof(TableauLobattoIIIAIIIB(2))     <: TableauPRK
    @test typeof(TableauLobattoIIIBIIIA(2))     <: TableauPRK

    # test instatiation of partioned tableau by composition of two RK tableaus
    @test typeof(TableauPRK(:PERK4, TableauRK4(), TableauRK4())) <: TableauPRK
    @test TableauPRK(:PERK4, TableauRK4(), TableauRK4()) == TableauPRK(:PERK4, TableauRK4())


    # TODO Add tests for TableauSARK, TableauSPARK and TableauGLM.

end
