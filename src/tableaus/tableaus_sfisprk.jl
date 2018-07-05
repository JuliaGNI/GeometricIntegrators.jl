# Tableau for the 2-stage stochastic LobattoIIIA-IIIB-IIID method (due to L. Jay)
function getTableauStochasticLobIIIABD2()

    TableauSFISPRK(:StochasticLobIIIABD2,getCoefficientsLobIIIA2(),getCoefficientsLobIIIA2(),
                                         getCoefficientsLobIIIB2(),getCoefficientsLobIIID2(),
                                         getCoefficientsLobIIIB2(),getCoefficientsLobIIID2())
end
