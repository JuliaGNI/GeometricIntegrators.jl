# Tableau for the 2-stage stochastic LobattoIIA-IIB method (Stormer-Verlet)
function getTableauStochasticStormerVerlet()

    TableauSFIPRK(:StochasticStormerVerlet,getCoefficientsLobIIIA2(),getCoefficientsLobIIIA2(),getCoefficientsLobIIIB2(),getCoefficientsLobIIIB2())
end
