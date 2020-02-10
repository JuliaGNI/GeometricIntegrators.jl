
@testset "$(rpad("Lagrange Polynomials Tests",80))" begin

    b = LagrangeBasis([0.0, 1.0])
    l = LagrangePolynomial([0.0, 1.0], ones(2))

    @test l.b == b
    @test l.c == ones(2)

    xp = [0.0, 0.5, 1.0]
    yp = zeros(3)
    yr = [1.0, 0.5, 0.0]

    evaluate!(l, xp, yp)
    @test yp == ones(3)

    y = zeros(2)
    y[1] = 1

    evaluate!(similar(l, y), xp, yp)
    @test yp == yr

    evaluate!(b, y, xp, yp)
    @test yp == yr

    y = zeros(2)
    y[2] = 1

    evaluate!(similar(l, y), xp, yp)
    @test yp == xp

    evaluate!(b, y, xp, yp)
    @test yp == xp

end
