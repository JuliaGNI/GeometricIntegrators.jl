
"Tableau for variational Gauss-Lobatto IIIA-IIIB method with two stages"
function getTableauLobIIIAIIIB2()
    d = [+1.0, -1.0]
    R∞ = -1
    TableauVPRK(:LobIIIAIIIB2, 2, getCoefficientsLobIIIA2(), getCoefficientsLobIIIB2(), R∞, d)
end

"Tableau for variational Gauss-Lobatto IIIA-IIIB method with three stages"
function getTableauLobIIIAIIIB3()
    d = [+1.0, -2.0, +1.0]
    R∞ = +1
    TableauVPRK(:LobIIIAIIIB3, 4, getCoefficientsLobIIIA3(), getCoefficientsLobIIIB3(), R∞, d)
end

"Tableau for variational Gauss-Lobatto IIIA-IIIB method with four stages"
function getTableauLobIIIAIIIB4()
    d = [+1.0, -√5, +√5, -1.0]
    R∞ = -1
    TableauVPRK(:LobIIIAIIIB4, 6, getCoefficientsLobIIIA4(), getCoefficientsLobIIIB4(), R∞, d)
end


"Tableau for variational Gauss-Lobatto IIIC-III method with two stages"
function getTableauLobIIIC2()
    R∞ = -1
    TableauVPRK(:LobIIIC2, 2, getCoefficientsLobIIIC2(), getCoefficientsLobIII2(), R∞)
end

"Tableau for variational Gauss-Lobatto IIIC-III method with three stages"
function getTableauLobIIIC3()
    R∞ = +1
    TableauVPRK(:LobIIIC3, 4, getCoefficientsLobIIIC3(), getCoefficientsLobIII3(), R∞)
end

"Tableau for variational Gauss-Lobatto IIIC-III method with four stages"
function getTableauLobIIIC4()
    R∞ = -1
    TableauVPRK(:LobIIIC4, 6, getCoefficientsLobIIIC4(), getCoefficientsLobIII4(), R∞)
end


"Tableau for variational Gauss-Lobatto IIID method with two stages"
function getTableauLobIIID2()
    lobD = getCoefficientsLobIIID2()
    R∞ = +1
    TableauVPRK(:LobIIID2, 2, lobD, lobD, R∞)
end

"Tableau for variational Gauss-Lobatto IIID method with three stages"
function getTableauLobIIID3()
    lobD = getCoefficientsLobIIID3()
    R∞ = -1
    TableauVPRK(:LobIIID3, 4, lobD, lobD, R∞)
end

"Tableau for variational Gauss-Lobatto IIID method with four stages"
function getTableauLobIIID4()
    lobD = getCoefficientsLobIIID4()
    R∞ = +1
    TableauVPRK(:LobIIID4, 6, lobD, lobD, R∞)
end


"Tableau for variational Gauss-Lobatto IIIE method with two stages"
function getTableauLobIIIE2()
    lobE = getCoefficientsLobIIIE2()
    R∞ = +1
    TableauVPRK(:LobIIIE2, 2, lobE, lobE, R∞)
end

"Tableau for variational Gauss-Lobatto IIIE method with three stages"
function getTableauLobIIIE3()
    lobE = getCoefficientsLobIIIE3()
    R∞ = -1
    TableauVPRK(:LobIIIE3, 4, lobE, lobE, R∞)
end

"Tableau for variational Gauss-Lobatto IIIE method with four stages"
function getTableauLobIIIE4()
    lobE = getCoefficientsLobIIIE4()
    R∞ = +1
    TableauVPRK(:LobIIIE4, 6, lobE, lobE, R∞)
end


"Tableau for variational Gauss-Lobatto IIIF method with two stages"
function getTableauLobIIIF2()
    lobF = getCoefficientsLobIIIF2()
    R∞ = +1
    TableauVPRK(:LobIIIF2, lobF.o, lobF, get_symplectic_conjugate_coefficients(lobF), R∞)
end

"Tableau for variational Gauss-Lobatto IIIF method with three stages"
function getTableauLobIIIF3()
    lobF = getCoefficientsLobIIIF3()
    R∞ = -1
    TableauVPRK(:LobIIIF3, lobF.o, lobF, get_symplectic_conjugate_coefficients(lobF), R∞)
end

"Tableau for variational Gauss-Lobatto IIIF method with four stages"
function getTableauLobIIIF4()
    lobF = getCoefficientsLobIIIF4()
    R∞ = +1
    TableauVPRK(:LobIIIF4, lobF.o, lobF, get_symplectic_conjugate_coefficients(lobF), R∞)
end


"Tableau for variational Gauss-Legendre method with s stages"
function getTableauVPGLRK(s; T=Float64)
    glrk = getCoefficientsGLRK(s, T=T)
    R∞ = (-1)^s
    TableauVPRK(Symbol("vpglrk", s), glrk.o, glrk, glrk, R∞)
end

"Tableau for variational symmetric Runge-Kutta method with 3 stages"
function getTableauVPSRK3()
    srk = getCoefficientsSRK3()
    R∞ = (-1)^3
    TableauVPRK(:vpsrk3, srk.o, srk, srk, R∞)
end
