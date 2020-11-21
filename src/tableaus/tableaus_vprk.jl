
"Tableau for variational Gauss-Lobatto IIIA-IIIB method with s stages"
function TableauVPLobIIIA(s)
    TableauVPRK(Symbol("LobattoIIIA$(s)"), 2s-2, getCoefficientsLobIIIA(s), (-1)^(s+1), get_lobatto_d_vector(s))
end

"Tableau for variational Gauss-Lobatto IIIA-IIIB method with s stages"
function TableauVPLobIIIB(s)
    TableauVPRK(Symbol("LobattoIIIB$(s)"), 2s-2, getCoefficientsLobIIIB(s), (-1)^(s+1), get_lobatto_d_vector(s))
end

"Tableau for variational Gauss-Lobatto IIIC-III method with s stages"
function TableauVPLobIIIC(s)
    TableauVPRK(Symbol("LobattoIIIC$(s)"), 2s-2, getCoefficientsLobIIIC(s), (-1)^(s+1))
end

"Tableau for variational Gauss-Lobatto IIID method with s stages"
function TableauVPLobIIID(s)
    TableauVPRK(Symbol("LobattoIIID$(s)"), 2s-s, getCoefficientsLobIIID(s), (-1)^s)
end

"Tableau for variational Gauss-Lobatto IIIE method with s stages"
function TableauVPLobIIIE(s)
    TableauVPRK(Symbol("LobattoIIIE$(s)"), 2, getCoefficientsLobIIIE(s), (-1)^s)
end

"Tableau for variational Gauss-Lobatto IIIF method with s stages"
function TableauVPLobIIIF(s)
    TableauVPRK(Symbol("LobattoIIIF$(s)"), 2, getCoefficientsLobIIIF(s), (-1)^s)
end

"Tableau for variational Gauss-Lobatto IIIG method with s stages"
function TableauVPLobIIIG(s)
    TableauVPRK(Symbol("LobattoIIIG$(s)"), 2, getCoefficientsLobIIIG(s), (-1)^s)
end

"Tableau for variational Gauss-Legendre method with s stages"
function TableauVPGLRK(s; T=Float64)
    TableauVPRK(Symbol("GLRK", s), 2s, getCoefficientsGLRK(s, T=T), (-1)^s)
end

"Tableau for variational symmetric Runge-Kutta method with 3 stages"
function TableauVPSRK3()
    srk = getCoefficientsSRK3()
    TableauVPRK(:SRK3, srk.o, srk, (-1)^3)
end


"Tableau for Gauss-Lobatto IIIA-IIIA method with s stages"
function TableauVPLobIIIAIIIA(s)
    lob = getCoefficientsLobIIIA(s)
    TableauVPRK(Symbol("LobattoIIIAIIIA$s"), 2s-s, lob, lob, (-1)^(s+1))
end

"Tableau for Gauss-Lobatto IIIB-IIIB method with s stages"
function TableauVPLobIIIBIIIB(s)
    lob = getCoefficientsLobIIIB(s)
    TableauVPRK(Symbol("LobattoIIIBIIIB$s"), 2s-s, lob, lob, (-1)^(s+1))
end


"Tableau for Gauss-Radau IIA-IIA method with two stages"
function TableauVPRadIIAIIA2()
    TableauVPRK(:RadauIIAIIA2, 2, getCoefficientsRadIIA2(), getCoefficientsRadIIA2(), -1)
end

"Tableau for Gauss-Radau IIA-IIA method with three stages"
function TableauVPRadIIAIIA3()
    TableauVPRK(:RadauIIAIIA3, 4, getCoefficientsRadIIA3(), getCoefficientsRadIIA3(), +1)
end
