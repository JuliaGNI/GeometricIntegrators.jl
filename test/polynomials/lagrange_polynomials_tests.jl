
lag = LagrangePolynomial([0.0, 1.0], ones(2))

xp = [0.0, 0.5, 1.0]
yp = zeros(3)
yr = [1.0, 0.5, 0.0]

y = zeros(2)
y[1] = 1
evaluate!(lag, y, xp, yp)
@test yp == yr

y = zeros(2)
y[2] = 1
evaluate!(lag, y, xp, yp)
@test yp == xp
