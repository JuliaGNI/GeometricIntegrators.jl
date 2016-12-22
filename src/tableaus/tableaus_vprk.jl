
function getTableauLobIIIA()
    a = [[0.0 0.0]
         [0.5 0.5]]
    b = [0.5, 0.5]
    c = [0.0, 1.0]

    TableauRK(:LobIIIA, 2, a, b, c)
end

function getTableauLobIIIB()
    a = [[0.5 0.0]
         [0.5 0.0]]
    b = [0.5, 0.5]
    c = [0.0, 1.0]

    TableauRK(:LobIIIB, 2, a, b, c)
end

"Tableau for variational Gauss-Lobatto IIIA-IIIB method with two stages"
function getTableauLobIIIAB2()
    d = [+1.0, -1.0]

    TableauVPRK(:LobIIIAB2, 2, getTableauLobIIIA(), getTableauLobIIIB(), d)
end
