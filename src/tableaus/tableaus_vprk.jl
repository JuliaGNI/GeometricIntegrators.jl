
"Tableau for variational Gauss-Lobatto IIIA-IIIB method with two stages"
function getTableauLobIIIAB2()
    a_q = [[0.0 0.0]
           [0.5 0.5]]
    b_q = [0.5, 0.5]
    c_q = [0.0, 1.0]

    a_p = [[0.5 0.0]
           [0.5 0.0]]
    b_p = [0.5, 0.5]
    c_p = [0.0, 1.0]

    d = [+1.0, -1.0]

    o = 2

    TableauVPRK(:LobIIIAB2, o, a_q, a_p, b_q, b_p, c_q, c_p, d)
end
