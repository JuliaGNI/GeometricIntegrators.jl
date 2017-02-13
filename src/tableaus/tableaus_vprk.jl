
"Tableau for variational Gauss-Lobatto IIIA-IIIB method with two stages"
function getTableauLobIIIAIIIB2()
    d = [+1.0, -1.0]
    R∞ = -1

    TableauVPRK(:LobIIIAB2, 2, getCoefficientsLobIIIA2(), getCoefficientsLobIIIB2(), R∞, d)
end

"Tableau for variational Gauss-Lobatto IIIA-IIIB method with three stages"
function getTableauLobIIIAIIIB3()
    d = [+0.5, -1.0, +0.5]
    R∞ = +1

    TableauVPRK(:LobIIIAB3, 4, getCoefficientsLobIIIA3(), getCoefficientsLobIIIB3(), R∞, d)
end

"Tableau for variational Gauss-Legendre method with s stages"
function getTableauVPGLRK(s)
    glrk = getCoefficientsGLRK(s)
    R∞ = -1^s

    TableauVPRK(Symbol("vpglrk", s), glrk.o, glrk, glrk, R∞)
end
