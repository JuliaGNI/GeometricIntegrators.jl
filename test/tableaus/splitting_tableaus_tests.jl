
@testset "$(rpad("Splitting Tableaus",80))" begin

    # instatiate all splitting tableaus
    @test typeof(TableauStrangA())       <: TableauSplitting{Float64}
    @test typeof(TableauStrangB())       <: TableauSplitting{Float64}

    @test typeof(TableauLieA())          <: TableauSplittingNS{Float64}
    @test typeof(TableauLieB())          <: TableauSplittingNS{Float64}
    @test typeof(TableauStrang())        <: TableauSplittingNS{Float64}
    @test typeof(TableauMcLachlan2())    <: TableauSplittingNS{Float64}
    @test typeof(TableauMcLachlan4())    <: TableauSplittingNS{Float64}

    @test typeof(TableauTripleJump())    <: TableauSplittingSS{Float64}
    @test typeof(TableauSuzukiFractal()) <: TableauSplittingSS{Float64}

end
