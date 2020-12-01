
"Implicit Euler"
function TableauImplicitEuler()
    a = ones(Float64, 1, 1)
    b = [1.0]
    c = [1.0]
    o = 1

    TableauFIRK(:implicit_euler, o, a, b, c)
end

"Implicit Midpoint"
function TableauImplicitMidpoint()
    a = 0.5*ones(Float64, 1, 1)
    b = [1.0]
    c = [0.5]
    o = 2

    TableauFIRK(:implicit_midpoint, o, a, b, c)
end

"Gauss-Legendre Runge-Kutta"
function TableauGLRK(s::Int)
    TableauFIRK(CoefficientsGLRK(s))
end

"Gauss-Lobatto-IIIA Runge-Kutta tableau with s stages."
function TableauLobattoIIIA(s)
    TableauFIRK(CoefficientsLobattoIIIA(s))
end

"Gauss-Lobatto-IIIB Runge-Kutta tableau with s stages."
function TableauLobattoIIIB(s)
    TableauFIRK(CoefficientsLobattoIIIB(s))
end

"Gauss-Lobatto-IIIC Runge-Kutta tableau with s stages."
function TableauLobattoIIIC(s)
    TableauFIRK(CoefficientsLobattoIIIC(s))
end
"Gauss-Lobatto-IIIC̄ Runge-Kutta tableau with s stages."
function TableauLobattoIIIC̄(s)
    TableauFIRK(CoefficientsLobattoIIIC̄(s))
end

"Gauss-Lobatto-IIID Runge-Kutta tableau with s stages."
function TableauLobattoIIID(s)
    TableauFIRK(CoefficientsLobattoIIID(s))
end

"Gauss-Lobatto-IIIE Runge-Kutta tableau with s stages."
function TableauLobattoIIIE(s)
    TableauFIRK(CoefficientsLobattoIIIE(s))
end

"Gauss-Lobatto-IIIF Runge-Kutta tableau with s stages."
function TableauLobattoIIIF(s)
    TableauFIRK(CoefficientsLobattoIIIF(s))
end

"Gauss-Lobatto-IIIG Runge-Kutta tableau with s stages."
function TableauLobattoIIIG(s)
    TableauFIRK(CoefficientsLobattoIIIG(s))
end

"Gauss-Radau-IIA Runge-Kutta tableau with s=2 stages."
function TableauRadauIA(s)
    TableauFIRK(CoefficientsRadauIA(s))
end

"Gauss-Radau-IIA Runge-Kutta tableau with s=3 stages."
function TableauRadauIIA(s)
    TableauFIRK(CoefficientsRadauIIA(s))
end

"Symmetric Runge-Kutta tableau with three stages."
function TableauSRK3()
    TableauFIRK(CoefficientsSRK3())
end
