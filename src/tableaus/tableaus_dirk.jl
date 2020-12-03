

"Tableau of Crank-Nicolson two-stage, 2nd order method"
function TableauCrankNicolson()
    a = [[ 0.0   0.0 ]
         [ 0.5   0.5 ]]
    b =  [ 0.5,  0.5 ]
    c =  [ 0.0,  1.0 ]
    o = 2

    TableauDIRK(:cranknicolson, o, a, b, c)
end

"Tableau of Kraaijevanger and Spijker's two-stage, 2nd order method"
function TableauKraaijevangerSpijker()
    a = [[ 0.5   0.0  ]
         [-0.5   2.0 ]]
    b =  [-0.5,  1.5 ]
    c =  [ 0.5,  1.5 ]
    o = 2

    TableauDIRK(:kraaijevangerspijker, o, a, b, c)
end

"Tableau of Qin and Zhang's symplectic two-stage, 2nd order method"
function TableauQinZhang()
    a = [[ 0.25   0.00 ]
         [ 0.50   0.25 ]]
    b =  [ 0.50,  0.50 ]
    c =  [ 0.25,  0.75 ]
    o = 2

    TableauDIRK(:qinzhang, o, a, b, c)
end

"Tableau of Crouzeix's two-stage, 3rd order method"
function TableauCrouzeix()
    fac = 0.5/√3
    a = [[ 0.5+fac 0.0    ]
         [-2.0*fac 0.5+fac]]
    b =  [0.5,     0.5    ]
    c =  [0.5+fac, 0.5-fac]
    o = 3

    TableauDIRK(:crouzeix, o, a, b, c)
end
