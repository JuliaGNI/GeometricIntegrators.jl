
"Implicit Midpoint"
function getTableauImplicitMidpoint()
    a = 0.5*ones(Float64, 1, 1)
    b = [1.0]
    c = [0.5]
    o = 2

    TableauNLIRK(:implicit_midpoint, o, a, b, c)
end

"Gauss-Legendre Runge-Kutta, s=1"
function getTableauGLRK1()
    a = 0.5*ones(Float64, 1, 1)
    b = [1.0]
    c = [0.5]
    o = 2

    TableauNLIRK(:glrk1, o, a, b, c)
end

"Gauss-Legendre Runge-Kutta, s=2"
function getTableauGLRK2()
    fac = √3/6.
    a = [[0.25     0.25-fac]
         [0.25+fac 0.25    ]]
    b = [0.5,     0.5    ]
    c = [0.5-fac, 0.5+fac]
    o = 4

    TableauNLIRK(:glrk2, o, a, b, c)
end
