"Tableau for the s-stage Gauss-Lobatto SFIRK method"
function getTableauStochasticGLRK(s::Int)

    if s==1
        name = :StochasticGLRK1
    elseif s==2
        name = :StochasticGLRK2
    elseif s==3
        name = :StochasticGLRK3
    elseif s==4
        name = :StochasticGLRK4
    elseif s==5
        name = :StochasticGLRK5
    else s==6
        name = :StochasticGLRK6
    end

    TableauSFIRK(name,getCoefficientsGLRK(s),getCoefficientsGLRK(s))
end


"""
Tableau for the 2-stage stochastic symplectic DIRK method
  Tableau for the stochastic symplectic DIRK method
  Satisfies the conditions for Lagrange-d'Alembert integrators.
  Satisfies the conditions for strong convergence of order 1.0 for one Wiener process
"""
function getTableauStochasticDIRK(c::Number=0.5)

    a_drift = [[c/2  0.0    ]
                [c    (1-c)/2]]

    b_drift = [c, 1-c]

    c_drift = [c/2, (1+c)/2]

    TableauSFIRK(:StochasticSymplecticDIRK, 2, a_drift, b_drift, c_drift, 2, a_drift, b_drift, c_drift)
end
