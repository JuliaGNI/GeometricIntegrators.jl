
@testset "$(rpad("Stochastic Runge-Kutta Tableaus",80))" begin

    @test typeof(getTableauStochasticEuler()) <: TableauSERK
    @test typeof(getTableauStochasticHeun()) <: TableauSERK
    @test typeof(getTableauPlaten()) <: TableauSERK
    @test typeof(getTableauBurrageR2()) <: TableauSERK
    @test typeof(getTableauBurrageCL()) <: TableauSERK
    @test typeof(getTableauBurrageE1()) <: TableauSERK
    @test typeof(getTableauBurrageG5()) <: TableauSERK

    @test typeof(getTableauStochasticGLRK(1)) <: TableauSIRK
    @test typeof(getTableauStochasticGLRK(2)) <: TableauSIRK
    @test typeof(getTableauStochasticDIRK()) <: TableauSIRK

    @test typeof(getTableauStochasticSymplecticEuler()) <: TableauSIPRK
    @test typeof(getTableauStochasticStormerVerlet()) <: TableauSIPRK

    @test typeof(getTableauStochasticLobIIIABD2()) <: TableauSISPRK
    @test typeof(getTableauModifiedStochasticStormerVerlet()) <: TableauSISPRK

    @test typeof(getTableauRosslerRS1()) <: TableauWERK
    @test typeof(getTableauRosslerRS2()) <: TableauWERK

    @test typeof(getTableauSRKw1()) <: TableauWIRK
    @test typeof(getTableauSRKw2()) <: TableauWIRK

end
