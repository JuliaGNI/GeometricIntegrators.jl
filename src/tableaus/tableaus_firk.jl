
"Implicit Euler"
function getTableauImplicitEuler()
    a = ones(Float64, 1, 1)
    b = [1.0]
    c = [1.0]
    o = 1

    TableauFIRK(:implicit_euler, o, a, b, c)
end

"Implicit Midpoint"
function getTableauImplicitMidpoint()
    a = 0.5*ones(Float64, 1, 1)
    b = [1.0]
    c = [0.5]
    o = 2

    TableauFIRK(:implicit_midpoint, o, a, b, c)
end

"Gauss-Legendre Runge-Kutta, s=1"
function getTableauGLRK1()
    a = 0.5*ones(Float64, 1, 1)
    b = [1.0]
    c = [0.5]
    o = 2

    TableauFIRK(:glrk1, o, a, b, c)
end

"Gauss-Legendre Runge-Kutta, s=2"
function getTableauGLRK2()
    a = [[1/4       1/4-√3/6]
         [1/4+√3/6  1/4     ]]
    b =  [1/2,      1/2     ]
    c =  [1/2-√3/6, 1/2+√3/6]
    o = 4

    TableauFIRK(:glrk2, o, a, b, c)
end

"Gauss-Legendre Runge-Kutta, s=3"
function getTableauGLRK3()
    a = [[5/36         2/9-√15/15  5/36-√15/30]
         [5/36+√15/24  2/9         5/36-√15/24]
         [5/36+√15/30  2/9+√15/15  5/36       ]]
    b =  [5/18,        4/9,        5/18       ]
    c =  [1/2-√15/10,  1/2,        1/2+√15/10 ]
    o = 6

    TableauFIRK(:glrk3, o, a, b, c)
end

"Gauss-Legendre Runge-Kutta, s=3"
function getTableauSRK3()
    TableauFIRK(getCoefficientsSRK3())
end

function getTableauGLRK(s::Int)
    TableauFIRK(getCoefficientsGLRK(s))
end
