
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

"Tableau for variational Gauss-Legendre method with s stages"
function getTableauVPGLRK(s)
    glrk = getCoefficientsGLRK(s)
    R∞ = (-1)^s

    TableauVPRK(Symbol("vpglrk", s), glrk.o, glrk, glrk, R∞)
end

"Tableau for variational symmetric Runge-Kutta method with 3 stages"
function getTableauVPSRK3()
    srk = getCoefficientsSRK3()
    R∞ = (-1)^3

    TableauVPRK(:vpsrk3, srk.o, srk, srk, R∞)
end
