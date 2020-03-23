
@testset "$(rpad("Lagrange Basis Tests",80))" begin

    x = [0.0, 1.0]
    b = LagrangeBasis(x)

    @test nodes(b) == x

    @test evaluate(b, 1, 0.0) == 1.0
    @test evaluate(b, 1, 0.5) == 0.5
    @test evaluate(b, 1, 1.0) == 0.0

    @test evaluate(b, 2, 0.0) == 0.0
    @test evaluate(b, 2, 0.5) == 0.5
    @test evaluate(b, 2, 1.0) == 1.0

    @test derivative(b, 1, 0.0) == -1.0
    @test derivative(b, 1, 0.5) == -1.0
    @test derivative(b, 1, 1.0) == -1.0

    @test derivative(b, 2, 0.0) == 1.0
    @test derivative(b, 2, 0.5) == 1.0
    @test derivative(b, 2, 1.0) == 1.0

    @test derivative(b, 1, 1) == -1.0
    @test derivative(b, 1, 2) == -1.0

    @test derivative(b, 2, 1) == 1.0
    @test derivative(b, 2, 2) == 1.0

    @test integral(b, 1, 0.0) == 0.0
    @test integral(b, 1, 0.5) == 0.375
    @test integral(b, 1, 1.0) == 0.5

    @test integral(b, 2, 0.0) == 0.0
    @test integral(b, 2, 0.5) == 0.125
    @test integral(b, 2, 1.0) == 0.5

    @test integral(b, 1, 1) == 0.0
    @test integral(b, 1, 2) == 0.5

    @test integral(b, 2, 1) == 0.0
    @test integral(b, 2, 2) == 0.5

end
