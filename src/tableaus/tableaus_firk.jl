
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
    TableauFIRK(getCoefficientsGLRK1())
end

"Gauss-Legendre Runge-Kutta, s=2"
function getTableauGLRK2()
    TableauFIRK(getCoefficientsGLRK2())
end

"Gauss-Legendre Runge-Kutta, s=3"
function getTableauGLRK3()
    TableauFIRK(getCoefficientsGLRK3())
end

"Gauss-Legendre Runge-Kutta, s=3"
function getTableauSRK3()
    TableauFIRK(getCoefficientsSRK3())
end

function getTableauGLRK(s::Int)
    TableauFIRK(getCoefficientsGLRK(s))
end
