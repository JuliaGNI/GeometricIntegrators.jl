
@doc raw"""
Tableau for variational Lobatto IIIA method with s stages

Uses Lobatto IIIA for the coefficients $a_{ij}$ and its symplectic conjugate, Lobatto IIIB, for the coefficients $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIA(::Type{T}, s) where {T}
    TableauVPRK(Symbol("LobattoIIIA$(s)"), 2s-2, TableauLobattoIIIA(T,s), (-1)^(s+1), get_lobatto_nullvector(s))
end

TableauVPLobattoIIIA(s) = TableauVPLobattoIIIA(Float64,s)


@doc raw"""
Tableau for variational Lobatto IIIB method with s stages

Uses Lobatto IIIB for the coefficients $a_{ij}$ and its symplectic conjugate, Lobatto IIIA, for the coefficients $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIB(::Type{T}, s) where {T}
    TableauVPRK(Symbol("LobattoIIIB$(s)"), 2s-2, TableauLobattoIIIB(T,s), (-1)^(s+1), get_lobatto_nullvector(s))
end

TableauVPLobattoIIIB(s) = TableauVPLobattoIIIB(Float64,s)


@doc raw"""
Tableau for variational Lobatto IIIC method with s stages

Uses Lobatto IIIC for the coefficients $a_{ij}$ and its symplectic conjugate, Lobatto IIIC̄, for the coefficients $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIC(::Type{T}, s) where {T}
    TableauVPRK(Symbol("LobattoIIIC$(s)"), 2s-2, TableauLobattoIIIC(T,s), (-1)^(s+1))
end

TableauVPLobattoIIIC(s) = TableauVPLobattoIIIC(Float64,s)


@doc raw"""
Tableau for variational Lobatto IIIC̄ method with s stages

Uses Lobatto IIIC̄ for the coefficients $a_{ij}$ and its symplectic conjugate, Lobatto IIIC, for the coefficients $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIC̄(::Type{T}, s) where {T}
    TableauVPRK(Symbol("LobattoIIIC̄$(s)"), 2s-2, TableauLobattoIIIC̄(T,s), (-1)^(s+1))
end

TableauVPLobattoIIIC̄(s) = TableauVPLobattoIIIC̄(Float64,s)


@doc raw"""
Tableau for variational Lobatto IIID method with s stages

Uses Lobatto IIID for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIID(::Type{T}, s) where {T}
    TableauVPRK(Symbol("LobattoIIID$(s)"), 2s-s, TableauLobattoIIID(T,s), (-1)^s)
end

TableauVPLobattoIIID(s) = TableauVPLobattoIIID(Float64,s)


@doc raw"""
Tableau for variational Lobatto IIIE method with s stages

Uses Lobatto IIIE for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIE(::Type{T}, s) where {T}
    TableauVPRK(Symbol("LobattoIIIE$(s)"), 2, TableauLobattoIIIE(T,s), (-1)^s)
end

TableauVPLobattoIIIE(s) = TableauVPLobattoIIIE(Float64,s)


@doc raw"""
Tableau for variational Lobatto IIIF method with s stages

Uses Lobatto IIIF for the coefficients $a_{ij}$ and its symplectic conjugate for the coefficients $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIF(::Type{T}, s) where {T}
    TableauVPRK(Symbol("LobattoIIIF$(s)"), 2, TableauLobattoIIIF(T,s), (-1)^s)
end

TableauVPLobattoIIIF(s) = TableauVPLobattoIIIF(Float64,s)


@doc raw"""
Tableau for variational Lobatto IIIG method with s stages

Uses Lobatto IIIG for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIG(::Type{T}, s) where {T}
    TableauVPRK(Symbol("LobattoIIIG$(s)"), 2, TableauLobattoIIIG(T,s), (-1)^s)
end

TableauVPLobattoIIIG(s) = TableauVPLobattoIIIG(Float64,s)


@doc raw"""
Tableau for variational Gauss-Legendre method with s stages

Uses Gauss coefficients for both $a_{ij}$ and $\bar{a}_{ij}$.
"""
function TableauVPGLRK(::Type{T}, s) where {T}
    TableauVPRK(Symbol("GLRK", s), 2s, TableauGauss(T, s), (-1)^s)
end

TableauVPGLRK(s) = TableauVPGLRK(Float64,s)


@doc raw"""
Tableau for variational symmetric Runge-Kutta method with 3 stages

Uses SRK3 coefficients for both $a_{ij}$ and $\bar{a}_{ij}$.
"""
function TableauVPSRK3(T=Float64)
    srk = TableauSRK3(T)
    TableauVPRK(:SRK3, srk.o, srk, (-1)^3)
end


@doc raw"""
Tableau for an implicit partitioned Lobatto IIIA Runge-Kutta method with s stages

Uses Lobatto IIIA for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIAIIIA(::Type{T}, s) where {T}
    lob = TableauLobattoIIIA(T,s)
    TableauVPRK(Symbol("LobattoIIIAIIIA$s"), 2s-2, lob, lob, (-1)^(s+1))
end

TableauVPLobattoIIIAIIIA(s) = TableauVPLobattoIIIAIIIA(Float64,s)


@doc raw"""
Tableau for an implicit partitioned Lobatto IIIB Runge-Kutta method with s stages

Uses Lobatto IIIB for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIBIIIB(::Type{T}, s) where {T}
    lob = TableauLobattoIIIB(T,s)
    TableauVPRK(Symbol("LobattoIIIBIIIB$s"), 2s-2, lob, lob, (-1)^(s+1))
end

TableauVPLobattoIIIBIIIB(s) = TableauVPLobattoIIIBIIIB(Float64,s)


@doc raw"""
Tableau for an implicit partitioned Radau IIA Runge-Kutta method with s stages

Uses Radau IIA for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
function TableauVPRadauIIAIIA(::Type{T}, s) where {T}
    rad = TableauRadauIIA(T,s)
    TableauVPRK(Symbol("RadauIIAIIA$s"), 2s-1, rad, rad, (-1)^(s+1))
end

TableauVPRadauIIAIIA(s) = TableauVPRadauIIAIIA(Float64,s)
