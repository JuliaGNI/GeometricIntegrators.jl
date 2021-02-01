
@doc raw"""
Tableau for variational Lobatto IIIA method with s stages

Uses Lobatto IIIA for the coefficients $a_{ij}$ and its symplectic conjugate, Lobatto IIIB, for the coefficients $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIA(s)
    TableauVPRK(Symbol("LobattoIIIA$(s)"), 2s-2, TableauLobattoIIIA(s), (-1)^(s+1), get_lobatto_d_vector(s))
end

@doc raw"""
Tableau for variational Lobatto IIIB method with s stages

Uses Lobatto IIIB for the coefficients $a_{ij}$ and its symplectic conjugate, Lobatto IIIA, for the coefficients $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIB(s)
    TableauVPRK(Symbol("LobattoIIIB$(s)"), 2s-2, TableauLobattoIIIB(s), (-1)^(s+1), get_lobatto_d_vector(s))
end

@doc raw"""
Tableau for variational Lobatto IIIC method with s stages

Uses Lobatto IIIC for the coefficients $a_{ij}$ and its symplectic conjugate, Lobatto IIIC̄, for the coefficients $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIC(s)
    TableauVPRK(Symbol("LobattoIIIC$(s)"), 2s-2, TableauLobattoIIIC(s), (-1)^(s+1))
end

@doc raw"""
Tableau for variational Lobatto IIIC̄ method with s stages

Uses Lobatto IIIC̄ for the coefficients $a_{ij}$ and its symplectic conjugate, Lobatto IIIC, for the coefficients $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIC̄(s)
    TableauVPRK(Symbol("LobattoIIIC̄$(s)"), 2s-2, CoefficientsLobattoIIIC̄(s), (-1)^(s+1))
end

@doc raw"""
Tableau for variational Lobatto IIID method with s stages

Uses Lobatto IIID for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIID(s)
    TableauVPRK(Symbol("LobattoIIID$(s)"), 2s-s, TableauLobattoIIID(s), (-1)^s)
end

@doc raw"""
Tableau for variational Lobatto IIIE method with s stages

Uses Lobatto IIIE for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIE(s)
    TableauVPRK(Symbol("LobattoIIIE$(s)"), 2, TableauLobattoIIIE(s), (-1)^s)
end

@doc raw"""
Tableau for variational Lobatto IIIF method with s stages

Uses Lobatto IIIF for the coefficients $a_{ij}$ and its symplectic conjugate for the coefficients $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIF(s)
    TableauVPRK(Symbol("LobattoIIIF$(s)"), 2, TableauLobattoIIIF(s), (-1)^s)
end

@doc raw"""
Tableau for variational Lobatto IIIG method with s stages

Uses Lobatto IIIG for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIG(s)
    TableauVPRK(Symbol("LobattoIIIG$(s)"), 2, TableauLobattoIIIG(s), (-1)^s)
end

@doc raw"""
Tableau for variational Gauss-Legendre method with s stages

Uses Gauss coefficients for both $a_{ij}$ and $\bar{a}_{ij}$.
"""
function TableauVPGLRK(s; T=Float64)
    TableauVPRK(Symbol("GLRK", s), 2s, TableauGauss(T, s), (-1)^s)
end

@doc raw"""
Tableau for variational symmetric Runge-Kutta method with 3 stages

Uses SRK3 coefficients for both $a_{ij}$ and $\bar{a}_{ij}$.
"""
function TableauVPSRK3()
    srk = TableauSRK3()
    TableauVPRK(:SRK3, srk.o, srk, (-1)^3)
end


@doc raw"""
Tableau for an implicit partitioned Lobatto IIIA Runge-Kutta method with s stages

Uses Lobatto IIIA for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIAIIIA(s)
    lob = TableauLobattoIIIA(s)
    TableauVPRK(Symbol("LobattoIIIAIIIA$s"), 2s-2, lob, lob, (-1)^(s+1))
end

@doc raw"""
Tableau for an implicit partitioned Lobatto IIIB Runge-Kutta method with s stages

Uses Lobatto IIIB for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIBIIIB(s)
    lob = TableauLobattoIIIB(s)
    TableauVPRK(Symbol("LobattoIIIBIIIB$s"), 2s-2, lob, lob, (-1)^(s+1))
end


@doc raw"""
Tableau for an implicit partitioned Radau IIA Runge-Kutta method with s stages

Uses Radau IIA for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
function TableauVPRadauIIAIIA(s)
    TableauVPRK(Symbol("RadauIIAIIA$s"), 2s-1, TableauRadauIIA(s), TableauRadauIIA(s), (-1)^(s+1))
end
