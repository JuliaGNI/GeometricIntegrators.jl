"""
Tableau for the 1-stage SRKw1 method due to Wang, Hong & Xu
  Method cited in
  Wang, Hong, Xu, "Construction of Symplectic Runge-Kutta Methods for Stochastic Hamiltonian Systems",
  Commun. Comput. Phys. 21(1), 2017
  According to the paper, the method has weak order 1.0.
"""
function TableauSRKw1(x::Number=0.0)

    A0 = 0.5*ones(typeof(x),1,1)

    A1 = (1. - x)*ones(typeof(x),1,1)

    B0 = x*ones(typeof(x),1,1)

    B1 = 0.5*ones(typeof(x),1,1)

    B3 = 0.5*ones(typeof(x),1,1)

    α  =  [1.0]
    β1 =  [1.0]

    c0 =  [0.5]
    c1 =  [1. - x]

    TableauWIRK(:SRKw1, A0, A1, B0, B1, B3, α, β1, c0, c1)
end


"""
Tableau for the 4-stage SRKw2 method due to Wang, Hong & Xu
  Method cited in
  Wang, Hong, Xu, "Construction of Symplectic Runge-Kutta Methods for Stochastic Hamiltonian Systems",
  Commun. Comput. Phys. 21(1), 2017
  According to the paper, the method has weak order 2.0 when applied to systems
  driven by one-dimensional noise.
"""
function TableauSRKw2(x1::Number=0.0, x2::Number=0.0, x3::Number=0.0)

    A0 = [[1. / 8.   0.0        0.0       0.0]
          [1. / 4.   1. / 8.    0.0       0.0]
          [1. / 4.   1. / 4.    1. / 8.   0.0]
          [1. / 4.   1. / 4.    1. / 4.   1. / 8.]]

    A1 = [[-1. / 6. + sqrt(3)/6.   1. / 3. - sqrt(3)/6.   0.0   1. / 3.]
          [ 1. / 2.                0.0                    0.0   0.0]
          [0.0                     0.0                    0.0   0.0]
          [0.0                     0.0                    0.0   0.0]]

    B0 = [[ 5. / 6. - sqrt(3)/3.     -1. / 2.    0.0   0.0]
          [-1. / 6. + sqrt(3)/3.      1. / 2.    0.0   0.0]
          [ 1. / 2.                   1. / 2.    0.0   0.0]
          [-1. / 6.                   1. / 2.    0.0   0.0]]

    B1 = [[1. / 4.                1. / 4. - sqrt(3)/6.   0.0   0.0]
          [1. / 4. + sqrt(3)/6.   1. / 4.                0.0   0.0]
          [x1                     x2                     0.0   x3]
          [0.0                    -1. / 2.               0.0   0.0]]

    B3 = zeros(typeof(x1),4,4)

    α  =  [0.25, 0.25, 0.25, 0.25]
    β1 =  [0.5,  0.5,  0.0,  0.0]

    c0 =  [1. / 8., 3. / 8., 5. / 8., 7. / 8.]
    c1 =  [0.5, 0.5, 0.0, 0.0]

    TableauWIRK(:SRKw2, A0, A1, B0, B1, B3, α, β1, c0, c1)
end
