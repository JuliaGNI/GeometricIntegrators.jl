"""
Tableau for the explicit 4-stage RS1 method due to Andreas Rossler
  Method cited in Table 5.2 in
  Andreas Rossler, "Second order Runge-Kutta methods for Stratonovich stochastic differential equations",
  BIT Numerical Mathematics (2007) 47
  According to the paper, the method has weak order 2.0.
"""
function getTableauRosslerRS1()

    A0 = [[0.0   0.0   0.0   0.0]
          [0.0   0.0   0.0   0.0]
          [1.    0.0   0.0   0.0]
          [0.0   0.0   0.0   0.0]]

    A1 = [[0.0   0.0   0.0   0.0]
          [0.0   0.0   0.0   0.0]
          [1.    0.0   0.0   0.0]
          [1.    0.0   0.0   0.0]]

    A2 =  zeros(Float64,4,4)

    B0 = [[0.0     0.0    0.0   0.0]
          [0.0     0.0    0.0   0.0]
          [1. / 4.   3. / 4.  0.0   0.0]
          [0.0     0.0    0.0   0.0]]

    B1 = [[0.0     0.0    0.0   0.0]
          [2. / 3.   0.0    0.0   0.0]
          [1. / 12.  1. / 4.  0.0   0.0]
          [-5. / 4.  1. / 4.  2.    0.0]]

    B2 = [[0.0     0.0    0.0   0.0]
          [1.      0.0    0.0   0.0]
          [-1.     0.0    0.0   0.0]
          [0.0     0.0    0.0   0.0]]

    B3 = [[0.0     0.0    0.0   0.0]
          [0.0     0.0    0.0   0.0]
          [1. / 4.   3. / 4.  0.0   0.0]
          [1. / 4.   3. / 4.  0.0   0.0]]

    α  =  [0.0, 0.0, 0.5, 0.5]
    β1 =  [1. / 8., 3. / 8., 3. / 8., 1. / 8.]
    β2 =  [0.0, -0.25, 0.25, 0.0]

    c0 =  [0.0, 0.0, 1., 0.0]
    c1 =  [0.0, 0.0, 1., 1.]
    c2 =  zeros(Float64, 4)

    TableauWERK(:RosslerRS1_explicit_method, A0, A1, A2, B0, B1, B2, B3, α, β1, β2, c0, c1, c2)
end


"""
Tableau for the explicit 4-stage RS2 method due to Andreas Rossler
  Method cited in Table 5.3 in
  Andreas Rossler, "Second order Runge-Kutta methods for Stratonovich stochastic differential equations",
  BIT Numerical Mathematics (2007) 47
  According to the paper, the method has weak order 2.0.
"""
function getTableauRosslerRS2()

    A0 = [[0.0      0.0      0.0   0.0]
          [2. / 3.  0.0      0.0   0.0]
          [1. / 6.  1. / 2.  0.0   0.0]
          [0.0      0.0      0.0   0.0]]

    A1 = [[0.0   0.0   0.0   0.0]
          [0.0   0.0   0.0   0.0]
          [1.    0.0   0.0   0.0]
          [1.    0.0   0.0   0.0]]

    A2 =  zeros(Float64,4,4)

    B0 = [[0.0       0.0      0.0   0.0]
          [0.0       0.0      0.0   0.0]
          [1. / 4.   3. / 4.  0.0   0.0]
          [0.0       0.0      0.0   0.0]]

    B1 = [[0.0       0.0      0.0   0.0]
          [2. / 3.   0.0      0.0   0.0]
          [1. / 12.  1. / 4.  0.0   0.0]
          [-5. / 4.  1. / 4.  2.    0.0]]

    B2 = [[0.0     0.0    0.0   0.0]
          [1.      0.0    0.0   0.0]
          [-1.     0.0    0.0   0.0]
          [0.0     0.0    0.0   0.0]]

    B3 = [[0.0       0.0      0.0   0.0]
          [0.0       0.0      0.0   0.0]
          [1. / 4.   3. / 4.  0.0   0.0]
          [1. / 4.   3. / 4.  0.0   0.0]]

    α  =  [0.25, 0.25, 0.5, 0.0]
    β1 =  [1. / 8., 3. / 8., 3. / 8., 1. / 8.]
    β2 =  [0.0, -0.25, 0.25, 0.0]

    c0 =  [0.0, 2. / 3., 2. / 3., 0.0]
    c1 =  [0.0, 0.0, 1., 1.]
    c2 =  zeros(Float64, 4)

    TableauWERK(:RosslerRS2_explicit_method, A0, A1, A2, B0, B1, B2, B3, α, β1, β2, c0, c1, c2)
end
