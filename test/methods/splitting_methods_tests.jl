
@testset "$(rpad("Splitting Methods",80))" begin

    # instatiate all splitting methods
    
    @test typeof(tableau(StrangA()))       <: TableauSplitting{Float64}
    @test typeof(tableau(StrangB()))       <: TableauSplitting{Float64}

    @test typeof(tableau(LieA()))          <: TableauSplittingNS{Float64}
    @test typeof(tableau(LieB()))          <: TableauSplittingNS{Float64}
    @test typeof(tableau(Strang()))        <: TableauSplittingNS{Float64}
    @test typeof(tableau(McLachlan2()))    <: TableauSplittingNS{Float64}
    @test typeof(tableau(McLachlan4()))    <: TableauSplittingNS{Float64}

    @test typeof(tableau(TripleJump()))    <: TableauSplittingSS{Float64}
    @test typeof(tableau(SuzukiFractal())) <: TableauSplittingSS{Float64}

end
