"Tableau for the s-stage Gauss-Lobatto SFIRK method"
function TableauStochasticGLRK(s::Int)
    TableauSIRK(Symbol("StochasticGLRK" * string(s)), CoefficientsGLRK(s), CoefficientsGLRK(s))
end


"""
Tableau for the 2-stage stochastic symplectic DIRK method
  Tableau for the stochastic symplectic DIRK method
  Satisfies the conditions for Lagrange-d'Alembert integrators.
  Satisfies the conditions for strong convergence of order 1.0 for one Wiener process
"""
function TableauStochasticDIRK(c::Number=0.5)

    a_drift = [[c/2  0.0    ]
                [c    (1-c)/2]]

    b_drift = [c, 1-c]

    c_drift = [c/2, (1+c)/2]

    TableauSIRK(:StochasticSymplecticDIRK, 2, a_drift, b_drift, c_drift, 2, a_drift, b_drift, c_drift)
end
