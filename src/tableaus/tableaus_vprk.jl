
"Tableau for variational Gauss-Lobatto IIIA-IIIB method with two stages"
function getTableauLobIIIAB2()
    d = [+1.0, -1.0]
    R∞ = -1

    TableauVPRK(:LobIIIAB2, 2, getCoefficientsLobIIIA(), getCoefficientsLobIIIB(), R∞, d)
end

"Tableau for variational Gauss-Legendre method with s stages"
function getTableauVPGLRK(s)
    glrk = getCoefficientsGLRK(s)
    R∞ = -1^s

    TableauVPRK(Symbol("vpglrk", s), glrk.o, glrk, glrk, R∞)
end
