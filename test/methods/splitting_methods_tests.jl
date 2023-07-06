
using GeometricIntegrators.Integrators: coefficients


@testset "$(rpad("Splitting Methods",80))" begin

    # instantiate all splitting methods
    
    @test typeof(coefficients(StrangA()))       <: SplittingCoefficientsGeneral{Float64}
    @test typeof(coefficients(StrangB()))       <: SplittingCoefficientsGeneral{Float64}

    @test typeof(coefficients(LieA()))          <: SplittingCoefficientsNonSymmetric{Float64}
    @test typeof(coefficients(LieB()))          <: SplittingCoefficientsNonSymmetric{Float64}
    @test typeof(coefficients(Strang()))        <: SplittingCoefficientsNonSymmetric{Float64}
    @test typeof(coefficients(McLachlan2()))    <: SplittingCoefficientsNonSymmetric{Float64}
    @test typeof(coefficients(McLachlan4()))    <: SplittingCoefficientsNonSymmetric{Float64}

    @test typeof(coefficients(TripleJump()))    <: SplittingCoefficientsSS{Float64}
    @test typeof(coefficients(SuzukiFractal())) <: SplittingCoefficientsSS{Float64}

end
