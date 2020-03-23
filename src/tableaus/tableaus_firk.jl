
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

"Gauss-Legendre Runge-Kutta"
function getTableauGLRK(s::Int)
    TableauFIRK(getCoefficientsGLRK(s))
end

"Gauss-Lobatto-IIIA Runge-Kutta, s=2"
function getTableauLobIIIA2()
    TableauFIRK(getCoefficientsLobIIIA2())
end

"Gauss-Lobatto-IIIA Runge-Kutta, s=3"
function getTableauLobIIIA3()
    TableauFIRK(getCoefficientsLobIIIA3())
end

"Gauss-Lobatto-IIIA Runge-Kutta, s=4"
function getTableauLobIIIA4()
    TableauFIRK(getCoefficientsLobIIIA4())
end

"Gauss-Lobatto-IIIB Runge-Kutta, s=2"
function getTableauLobIIIB2()
    TableauFIRK(getCoefficientsLobIIIB2())
end

"Gauss-Lobatto-IIIB Runge-Kutta, s=3"
function getTableauLobIIIB3()
    TableauFIRK(getCoefficientsLobIIIB3())
end

"Gauss-Lobatto-IIIB Runge-Kutta, s=4"
function getTableauLobIIIB4()
    TableauFIRK(getCoefficientsLobIIIB4())
end

"Gauss-Lobatto-IIIC Runge-Kutta, s=2"
function getTableauLobIIIC2()
    TableauFIRK(getCoefficientsLobIIIC2())
end

"Gauss-Lobatto-IIIC Runge-Kutta, s=3"
function getTableauLobIIIC3()
    TableauFIRK(getCoefficientsLobIIIC3())
end

"Gauss-Lobatto-IIIC Runge-Kutta, s=4"
function getTableauLobIIIC4()
    TableauFIRK(getCoefficientsLobIIIC4())
end

"Gauss-Lobatto-IIID Runge-Kutta, s=2"
function getTableauLobIIID2()
    TableauFIRK(getCoefficientsLobIIID2())
end

"Gauss-Lobatto-IIID Runge-Kutta, s=3"
function getTableauLobIIID3()
    TableauFIRK(getCoefficientsLobIIID3())
end

"Gauss-Lobatto-IIID Runge-Kutta, s=4"
function getTableauLobIIID4()
    TableauFIRK(getCoefficientsLobIIID4())
end

"Gauss-Lobatto-IIIE Runge-Kutta, s=2"
function getTableauLobIIIE2()
    TableauFIRK(getCoefficientsLobIIIE2())
end

"Gauss-Lobatto-IIIE Runge-Kutta, s=3"
function getTableauLobIIIE3()
    TableauFIRK(getCoefficientsLobIIIE3())
end

"Gauss-Lobatto-IIIE Runge-Kutta, s=4"
function getTableauLobIIIE4()
    TableauFIRK(getCoefficientsLobIIIE4())
end

"Gauss-Lobatto-IIIF Runge-Kutta, s=2"
function getTableauLobIIIF2()
    TableauFIRK(getCoefficientsLobIIIF2())
end

"Gauss-Lobatto-IIIF Runge-Kutta, s=3"
function getTableauLobIIIF3()
    TableauFIRK(getCoefficientsLobIIIF3())
end

"Gauss-Lobatto-IIIF Runge-Kutta, s=4"
function getTableauLobIIIF4()
    TableauFIRK(getCoefficientsLobIIIF4())
end

"Gauss-Lobatto-IIIG Runge-Kutta, s=2"
function getTableauLobIIIG2()
    TableauFIRK(getCoefficientsLobIIIG2())
end

"Gauss-Lobatto-IIIG Runge-Kutta, s=3"
function getTableauLobIIIG3()
    TableauFIRK(getCoefficientsLobIIIG3())
end

"Gauss-Lobatto-IIIG Runge-Kutta, s=4"
function getTableauLobIIIG4()
    TableauFIRK(getCoefficientsLobIIIG4())
end

"Gauss-Radau-IIA Runge-Kutta, s=2"
function getTableauRadIIA2()
    TableauFIRK(getCoefficientsRadIIA2())
end

"Gauss-Radau-IIA Runge-Kutta, s=3"
function getTableauRadIIA3()
    TableauFIRK(getCoefficientsRadIIA3())
end

"Gauss-Legendre Runge-Kutta, s=3"
function getTableauSRK3()
    TableauFIRK(getCoefficientsSRK3())
end
