
function getCoefficientsLobIIIA()
    a = [[0.0 0.0]
         [0.5 0.5]]
    b = [0.5, 0.5]
    c = [0.0, 1.0]

    CoefficientsRK(:LobIIIA, 2, a, b, c)
end

function getCoefficientsLobIIIB()
    a = [[0.5 0.0]
         [0.5 0.0]]
    b = [0.5, 0.5]
    c = [0.0, 1.0]

    CoefficientsRK(:LobIIIB, 2, a, b, c)
end

"Tableau for variational Gauss-Lobatto IIIA-IIIB method with two stages"
function getTableauLobIIIAB2()
    d = [+1.0, -1.0]
    R∞ = -1

    TableauVPRK(:LobIIIAB2, 2, getCoefficientsLobIIIA(), getCoefficientsLobIIIB(), R∞, d)
end

"Tableau for variational Gauss-Legendre method with s stages"
function getTableauVPGLRK(s)
    glrk = getTableauGLRK(s)
    R∞ = -1^s

    TableauVPRK(Symbol("vglrk", s), 2, glrk.q, glrk.q, R∞)
end
