
function TableauCrouzeix()
    fac = 0.5/√3
    a = [[ 0.5+fac 0.0    ]
         [-2.0*fac 0.5+fac]]
    b = [0.5,     0.5    ]
    c = [0.5+fac, 0.5-fac]
    o = 3

    TableauDIRK(:crouzeix, o, a, b, c)
end
