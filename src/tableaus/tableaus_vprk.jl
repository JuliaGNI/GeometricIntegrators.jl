
"Tableau for variational Gauss-Lobatto IIIA-IIIB method with two stages"
function getTableauVPLobIIIA2()
    R∞ = -1
    TableauVPRK(:LobIIIA2, 2, getCoefficientsLobIIIA2(), R∞, get_lobatto_d_vector(2))
end

"Tableau for variational Gauss-Lobatto IIIA-IIIB method with three stages"
function getTableauVPLobIIIA3()
    R∞ = +1
    TableauVPRK(:LobIIIA3, 4, getCoefficientsLobIIIA3(), R∞, get_lobatto_d_vector(3))
end

"Tableau for variational Gauss-Lobatto IIIA-IIIB method with four stages"
function getTableauVPLobIIIA4()
    R∞ = -1
    TableauVPRK(:LobIIIA4, 6, getCoefficientsLobIIIA4(), R∞, get_lobatto_d_vector(4))
end


"Tableau for variational Gauss-Lobatto IIIA-IIIB method with two stages"
function getTableauVPLobIIIB2()
    R∞ = -1
    TableauVPRK(:LobIIIB2, 2, getCoefficientsLobIIIB2(), R∞, get_lobatto_d_vector(2))
end

"Tableau for variational Gauss-Lobatto IIIA-IIIB method with three stages"
function getTableauVPLobIIIB3()
    R∞ = +1
    TableauVPRK(:LobIIIB3, 4, getCoefficientsLobIIIB3(), R∞, get_lobatto_d_vector(3))
end

"Tableau for variational Gauss-Lobatto IIIA-IIIB method with four stages"
function getTableauVPLobIIIB4()
    R∞ = -1
    TableauVPRK(:LobIIIB4, 6, getCoefficientsLobIIIB4(), R∞, get_lobatto_d_vector(4))
end


"Tableau for variational Gauss-Lobatto IIIC-III method with two stages"
function getTableauVPLobIIIC2()
    R∞ = -1
    TableauVPRK(:LobIIIC2, 2, getCoefficientsLobIIIC2(), R∞)
end

"Tableau for variational Gauss-Lobatto IIIC-III method with three stages"
function getTableauVPLobIIIC3()
    R∞ = +1
    TableauVPRK(:LobIIIC3, 4, getCoefficientsLobIIIC3(), R∞)
end

"Tableau for variational Gauss-Lobatto IIIC-III method with four stages"
function getTableauVPLobIIIC4()
    R∞ = -1
    TableauVPRK(:LobIIIC4, 6, getCoefficientsLobIIIC4(), R∞)
end


"Tableau for variational Gauss-Lobatto IIID method with two stages"
function getTableauVPLobIIID2()
    lobD = getCoefficientsLobIIID2()
    R∞ = +1
    TableauVPRK(:LobIIID2, 2, lobD, R∞)
end

"Tableau for variational Gauss-Lobatto IIID method with three stages"
function getTableauVPLobIIID3()
    lobD = getCoefficientsLobIIID3()
    R∞ = -1
    TableauVPRK(:LobIIID3, 4, lobD, R∞)
end

"Tableau for variational Gauss-Lobatto IIID method with four stages"
function getTableauVPLobIIID4()
    lobD = getCoefficientsLobIIID4()
    R∞ = +1
    TableauVPRK(:LobIIID4, 6, lobD, R∞)
end


"Tableau for variational Gauss-Lobatto IIIE method with two stages"
function getTableauVPLobIIIE2()
    lobE = getCoefficientsLobIIIE2()
    R∞ = +1
    TableauVPRK(:LobIIIE2, 2, lobE, R∞)
end

"Tableau for variational Gauss-Lobatto IIIE method with three stages"
function getTableauVPLobIIIE3()
    lobE = getCoefficientsLobIIIE3()
    R∞ = -1
    TableauVPRK(:LobIIIE3, 4, lobE, R∞)
end

"Tableau for variational Gauss-Lobatto IIIE method with four stages"
function getTableauVPLobIIIE4()
    lobE = getCoefficientsLobIIIE4()
    R∞ = +1
    TableauVPRK(:LobIIIE4, 6, lobE, R∞)
end


"Tableau for variational Gauss-Lobatto IIIF method with two stages"
function getTableauVPLobIIIF2()
    lobF = getCoefficientsLobIIIF2()
    R∞ = +1
    TableauVPRK(:LobIIIF2, lobF.o, lobF, R∞)
end

"Tableau for variational Gauss-Lobatto IIIF method with three stages"
function getTableauVPLobIIIF3()
    lobF = getCoefficientsLobIIIF3()
    R∞ = -1
    TableauVPRK(:LobIIIF3, lobF.o, lobF, R∞)
end

"Tableau for variational Gauss-Lobatto IIIF method with four stages"
function getTableauVPLobIIIF4()
    lobF = getCoefficientsLobIIIF4()
    R∞ = +1
    TableauVPRK(:LobIIIF4, lobF.o, lobF, R∞)
end


"Tableau for variational Gauss-Lobatto IIIG method with two stages"
function getTableauVPLobIIIG2()
    lobG = getCoefficientsLobIIIG2()
    R∞ = +1
    TableauVPRK(:LobIIIG2, 4, lobG, lobG, R∞)
end

"Tableau for variational Gauss-Lobatto IIIG method with three stages"
function getTableauVPLobIIIG3()
    lobG = getCoefficientsLobIIIG3()
    R∞ = -1
    TableauVPRK(:LobIIIG3, 6, lobG, lobG, R∞)
end

"Tableau for variational Gauss-Lobatto IIIG method with four stages"
function getTableauVPLobIIIG4()
    lobG = getCoefficientsLobIIIG4()
    R∞ = +1
    TableauVPRK(:LobIIIG4, 8, lobG, lobG, R∞)
end


"Tableau for variational Gauss-Legendre method with s stages"
function getTableauVPGLRK(s; T=Float64)
    glrk = getCoefficientsGLRK(s, T=T)
    R∞ = (-1)^s
    TableauVPRK(Symbol("GLRK", s), glrk.o, glrk, R∞)
end

"Tableau for variational symmetric Runge-Kutta method with 3 stages"
function getTableauVPSRK3()
    srk = getCoefficientsSRK3()
    R∞ = (-1)^3
    TableauVPRK(:SRK3, srk.o, srk, R∞)
end


"Tableau for Gauss-Lobatto IIIA-IIIA method with two stages"
function getTableauVPLobIIIAIIIA2()
    R∞ = -1
    TableauVPRK(:LobIIIAIIIA2, 2, getCoefficientsLobIIIA2(), getCoefficientsLobIIIA2(), R∞)
end

"Tableau for Gauss-Lobatto IIIA-IIIA method with three stages"
function getTableauVPLobIIIAIIIA3()
    R∞ = +1
    TableauVPRK(:LobIIIAIIIA3, 4, getCoefficientsLobIIIA3(), getCoefficientsLobIIIA3(), R∞)
end

"Tableau for Gauss-Lobatto IIIA-IIIA method with four stages"
function getTableauVPLobIIIAIIIA4()
    R∞ = -1
    TableauVPRK(:LobIIIAIIIA4, 6, getCoefficientsLobIIIA4(), getCoefficientsLobIIIA4(), R∞)
end



"Tableau for Gauss-Radau IIA-IIA method with two stages"
function getTableauVPRadIIAIIA2()
    R∞ = -1
    TableauVPRK(:RadIIAIIA2, 2, getCoefficientsRadIIA2(), getCoefficientsRadIIA2(), R∞)
end

"Tableau for Gauss-Radau IIA-IIA method with three stages"
function getTableauVPRadIIAIIA3()
    R∞ = +1
    TableauVPRK(:RadIIAIIA3, 4, getCoefficientsRadIIA3(), getCoefficientsRadIIA3(), R∞)
end
