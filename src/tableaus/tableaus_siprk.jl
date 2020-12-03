"""
Tableau for the stochastic symplectic Euler method
  Tableau for the stochastic symplectic Euler method
  Satisfies the conditions for Lagrange-d'Alembert integrators.
  Satisfies the conditions for strong convergence of order 1.0 for one Wiener process
  for special choices of the stochastic Hamiltonians and forces, e.g., h=h(q), f=0.
"""
function TableauStochasticSymplecticEuler()

    a_q = ones(Float64, 1, 1)
    b_q = [1.0]
    c_q = [1.0]

    a_p = zeros(Float64, 1, 1)
    b_p = [1.0]
    c_p = [0.0]

    TableauSIPRK(:StochasticSymplecticEuler, 1, a_q, b_q, c_q, 1, a_q, b_q, c_q, 1, a_p, b_p, c_p, 1, a_p, b_p, c_p)
end



"Tableau for the 2-stage stochastic LobattoIIA-IIB method (Stormer-Verlet)"
function TableauStochasticStormerVerlet()

    TableauSIPRK(:StochasticStormerVerlet,
                 CoefficientsLobattoIIIA(2), CoefficientsLobattoIIIA(2),
                 CoefficientsLobattoIIIB(2), CoefficientsLobattoIIIB(2))
end
