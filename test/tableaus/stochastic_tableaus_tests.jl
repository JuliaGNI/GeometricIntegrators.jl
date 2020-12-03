
@testset "$(rpad("Stochastic Runge-Kutta Tableaus",80))" begin

    @test typeof(TableauStochasticEuler()) <: TableauSERK
    @test typeof(TableauStochasticHeun()) <: TableauSERK
    @test typeof(TableauPlaten()) <: TableauSERK
    @test typeof(TableauBurrageR2()) <: TableauSERK
    @test typeof(TableauBurrageCL()) <: TableauSERK
    @test typeof(TableauBurrageE1()) <: TableauSERK
    @test typeof(TableauBurrageG5()) <: TableauSERK

    @test typeof(TableauStochasticGLRK(1)) <: TableauSIRK
    @test typeof(TableauStochasticGLRK(2)) <: TableauSIRK
    @test typeof(TableauStochasticDIRK()) <: TableauSIRK

    @test typeof(TableauStochasticSymplecticEuler()) <: TableauSIPRK
    @test typeof(TableauStochasticStormerVerlet()) <: TableauSIPRK

    @test typeof(TableauStochasticLobattoIIIABD2()) <: TableauSISPRK
    @test typeof(TableauModifiedStochasticStormerVerlet()) <: TableauSISPRK

    @test typeof(TableauRosslerRS1()) <: TableauWERK
    @test typeof(TableauRosslerRS2()) <: TableauWERK

    @test typeof(TableauSRKw1()) <: TableauWIRK
    @test typeof(TableauSRKw2()) <: TableauWIRK

end
