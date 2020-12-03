
"Tableau for variational Gauss-Lobatto IIIA-IIIB method with s stages"
function TableauVPLobattoIIIA(s)
    TableauVPRK(Symbol("LobattoIIIA$(s)"), 2s-2, CoefficientsLobattoIIIA(s), (-1)^(s+1), get_lobatto_d_vector(s))
end

"Tableau for variational Gauss-Lobatto IIIA-IIIB method with s stages"
function TableauVPLobattoIIIB(s)
    TableauVPRK(Symbol("LobattoIIIB$(s)"), 2s-2, CoefficientsLobattoIIIB(s), (-1)^(s+1), get_lobatto_d_vector(s))
end

"Tableau for variational Gauss-Lobatto IIIC-III method with s stages"
function TableauVPLobattoIIIC(s)
    TableauVPRK(Symbol("LobattoIIIC$(s)"), 2s-2, CoefficientsLobattoIIIC(s), (-1)^(s+1))
end

"Tableau for variational Gauss-Lobatto IIID method with s stages"
function TableauVPLobattoIIID(s)
    TableauVPRK(Symbol("LobattoIIID$(s)"), 2s-s, CoefficientsLobattoIIID(s), (-1)^s)
end

"Tableau for variational Gauss-Lobatto IIIE method with s stages"
function TableauVPLobattoIIIE(s)
    TableauVPRK(Symbol("LobattoIIIE$(s)"), 2, CoefficientsLobattoIIIE(s), (-1)^s)
end

"Tableau for variational Gauss-Lobatto IIIF method with s stages"
function TableauVPLobattoIIIF(s)
    TableauVPRK(Symbol("LobattoIIIF$(s)"), 2, CoefficientsLobattoIIIF(s), (-1)^s)
end

"Tableau for variational Gauss-Lobatto IIIG method with s stages"
function TableauVPLobattoIIIG(s)
    TableauVPRK(Symbol("LobattoIIIG$(s)"), 2, CoefficientsLobattoIIIG(s), (-1)^s)
end

"Tableau for variational Gauss-Legendre method with s stages"
function TableauVPGLRK(s; T=Float64)
    TableauVPRK(Symbol("GLRK", s), 2s, CoefficientsGLRK(s, T=T), (-1)^s)
end

"Tableau for variational symmetric Runge-Kutta method with 3 stages"
function TableauVPSRK3()
    srk = CoefficientsSRK3()
    TableauVPRK(:SRK3, srk.o, srk, (-1)^3)
end


"Tableau for Gauss-Lobatto IIIA-IIIA method with s stages"
function TableauVPLobattoIIIAIIIA(s)
    lob = CoefficientsLobattoIIIA(s)
    TableauVPRK(Symbol("LobattoIIIAIIIA$s"), 2s-2, lob, lob, (-1)^(s+1))
end

"Tableau for Gauss-Lobatto IIIB-IIIB method with s stages"
function TableauVPLobattoIIIBIIIB(s)
    lob = CoefficientsLobattoIIIB(s)
    TableauVPRK(Symbol("LobattoIIIBIIIB$s"), 2s-2, lob, lob, (-1)^(s+1))
end


"Tableau for Gauss-Radau IIA-IIA method with two stages"
function TableauVPRadauIIAIIA(s)
    TableauVPRK(Symbol("RadauIIAIIA$s"), 2s-1, CoefficientsRadauIIA(s), CoefficientsRadauIIA(s), (-1)^(s+1))
end
