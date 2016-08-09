
function getTableauExplicitMidpoint()
    a = [[0.0 0.0]
         [0.5 0.0]]
    b = [0.0, 1.0]
    c = [0.0, 0.5]
    o = 2

    TableauERK(:explicit_midpoint, o, a, b, c)
end

function getTableauHeun()
    a = [[0.0 0.0]
         [1.0 0.0]]
    b = [0.5, 0.5]
    c = [0.0, 1.0]
    o = 2

    TableauERK(:heun, o, a, b, c)
end
