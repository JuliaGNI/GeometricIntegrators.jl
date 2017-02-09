
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
