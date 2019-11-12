
@testset "$(rpad("Splitting Tableaus",80))" begin

    # instatiate all splitting tableaus
    @test typeof(getTableauLieA()) <: TableauSplittingNS{Float64}
    @test typeof(getTableauLieB()) <: TableauSplittingNS{Float64}
    @test typeof(getTableauStrang()) <: TableauSplittingNS{Float64}
    @test typeof(getTableauMcLachlan2()) <: TableauSplittingNS{Float64}
    @test typeof(getTableauMcLachlan4()) <: TableauSplittingNS{Float64}

    @test typeof(getTableauTripleJump()) <: TableauSplittingSS{Float64}
    @test typeof(getTableauSuzukiFractal()) <: TableauSplittingSS{Float64}

end
