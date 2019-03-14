"""
Tableau for the 2-stage stochastic LobattoIIIA-IIIB-IIID method
  Tableau for the 2-stage stochastic LobattoIIIA-IIIB-IIID method
  (based on the deterministic LobattoIIIA-IIIB-IIID due to L. Jay)
  It satisfies the conditions for convergence of order 1.0 for one Wiener process,
  but it doesn't satisfy the conditions for Lagrange-d'Alembert integrators
"""
function getTableauStochasticLobIIIABD2()

    TableauSFISPRK(:StochasticLobIIIABD2,getCoefficientsLobIIIA2(),getCoefficientsLobIIIA2(),
                                         getCoefficientsLobIIIB2(),getCoefficientsLobIIID2(),
                                         getCoefficientsLobIIIB2(),getCoefficientsLobIIID2())
end


"""
Tableau for the 2-stage modified stochastic LobattoIIIA-IIIB method
  Tableau for the 2-stage modified stochastic LobattoIIIA-IIIB method
  Satisfies the conditions for Lagrange-d'Alembert integrators
  and the conditions for convergence of order 1.0 for one Wiener process
"""
function getTableauModifiedStochasticStormerVerlet(c::Number=0.0)

    @assert c ≥ 0.0
    @assert c ≤ 1.0

    a_drift2 = [[c 0.0]
                [c 0.0]]

    b_drift2 = [c, 1-c]

    c_drift2 = [c, c]

    TableauSFISPRK(:StochasticModifiedStormerVerlet,getCoefficientsLobIIIA2(),getCoefficientsLobIIIA2(),
                                         getCoefficientsLobIIIB2(),CoefficientsRK(typeof(c), :cTableau, 1, a_drift2, b_drift2, c_drift2),
                                         getCoefficientsLobIIIB2(),getCoefficientsLobIIIB2())
end
