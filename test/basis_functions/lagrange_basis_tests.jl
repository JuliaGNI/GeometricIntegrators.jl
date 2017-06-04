
# check linear Lagrange basis
b = LagrangeBasis([0.0, 1.0])

@test lagrange(b, 1, 0.0) == 1.0
@test lagrange(b, 1, 0.5) == 0.5
@test lagrange(b, 1, 1.0) == 0.0

@test lagrange(b, 2, 0.0) == 0.0
@test lagrange(b, 2, 0.5) == 0.5
@test lagrange(b, 2, 1.0) == 1.0

@test lagrange_derivative(b, 1, 0.0) == -1.0
@test lagrange_derivative(b, 1, 0.5) == -1.0
@test lagrange_derivative(b, 1, 1.0) == -1.0

@test lagrange_derivative(b, 2, 0.0) == 1.0
@test lagrange_derivative(b, 2, 0.5) == 1.0
@test lagrange_derivative(b, 2, 1.0) == 1.0

@test lagrange_derivative(b, 1, 1) == -1.0
@test lagrange_derivative(b, 1, 2) == -1.0

@test lagrange_derivative(b, 2, 1) == 1.0
@test lagrange_derivative(b, 2, 2) == 1.0

@test lagrange_integral(b, 1, 0.0) == 0.0
@test lagrange_integral(b, 1, 0.5) == 0.375
@test lagrange_integral(b, 1, 1.0) == 0.5

@test lagrange_integral(b, 2, 0.0) == 0.0
@test lagrange_integral(b, 2, 0.5) == 0.125
@test lagrange_integral(b, 2, 1.0) == 0.5

@test lagrange_integral(b, 1, 1) == 0.0
@test lagrange_integral(b, 1, 2) == 0.5

@test lagrange_integral(b, 2, 1) == 0.0
@test lagrange_integral(b, 2, 2) == 0.5
