
@testset "$(rpad("Bernstein Basis Tests",80))" begin

    x = [0.0, 1.0]
    b = BernsteinBasis(x)

    @test nodes(b) == x

    @test evaluate(b, 1, 0.0) == 1.0
    @test evaluate(b, 1, 0.5) == 0.5
    @test evaluate(b, 1, 1.0) == 0.0

    @test evaluate(b, 2, 0.0) == 0.0
    @test evaluate(b, 2, 0.5) == 0.5
    @test evaluate(b, 2, 1.0) == 1.0

    @test deriv_basis(b, 1, 0.0) == -1.0
    @test deriv_basis(b, 1, 0.5) == -1.0
    @test deriv_basis(b, 1, 1.0) == -1.0

    @test deriv_basis(b, 2, 0.0) == 1.0
    @test deriv_basis(b, 2, 0.5) == 1.0
    @test deriv_basis(b, 2, 1.0) == 1.0

    @test deriv_basis(b, 1, 1) == -1.0
    @test deriv_basis(b, 1, 2) == -1.0

    @test deriv_basis(b, 2, 1) == 1.0
    @test deriv_basis(b, 2, 2) == 1.0

end
