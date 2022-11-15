
@doc raw"""
Tableau for variational Gauss-Legendre method with s stages

Uses Gauss coefficients for both $a_{ij}$ and $\bar{a}_{ij}$.
"""
function TableauVPGLRK(::Type{T}, s) where {T}
    tab = TableauGauss(T, s)
    TableauVPRK(Symbol("GLRK$s"), tab.o, tab, (-1)^s)
end

TableauVPGLRK(s) = TableauVPGLRK(Float64,s)


@doc raw"""
Tableau for variational symmetric Runge-Kutta method with 3 stages

Uses SRK3 coefficients for both $a_{ij}$ and $\bar{a}_{ij}$.
"""
function TableauVPSRK3(T=Float64)
    tab = TableauSRK3(T)
    TableauVPRK(:SRK3, tab.o, tab, (-1)^3)
end


@doc raw"""
Tableau for an implicit partitioned Radau IIA Runge-Kutta method with s stages

Uses Radau IIA for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
function TableauVPRadauIIA(::Type{T}, s) where {T}
    tab = TableauRadauIIA(T,s)
    TableauVPRK(Symbol("RadauIIA$s"), 2s-1, tab, tab, (-1)^(s+1))
end

TableauVPRadauIIA(s) = TableauVPRadauIIA(Float64,s)


@doc raw"""
Tableau for an implicit partitioned Radau IIB Runge-Kutta method with s stages

Uses Radau IIB for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
function TableauVPRadauIIB(::Type{T}, s) where {T}
    tab = TableauRadauIIB(T,s)
    TableauVPRK(Symbol("RadauIIB$s"), 2s-1, tab, tab, (-1)^(s+1))
end

TableauVPRadauIIB(s) = TableauVPRadauIIB(Float64,s)


@doc raw"""
Tableau for an implicit partitioned Lobatto IIIA Runge-Kutta method with s stages

Uses Lobatto IIIA for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIA(::Type{T}, s) where {T}
    lob = TableauLobattoIIIA(T,s)
    TableauVPRK(Symbol("LobattoIIIA$s"), 2s-2, lob, lob, (-1)^(s+1))
end

TableauVPLobattoIIIA(s) = TableauVPLobattoIIIA(Float64,s)


@doc raw"""
Tableau for an implicit partitioned Lobatto IIIB Runge-Kutta method with s stages

Uses Lobatto IIIB for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIB(::Type{T}, s) where {T}
    lob = TableauLobattoIIIB(T,s)
    TableauVPRK(Symbol("LobattoIIIB$s"), 2s-2, lob, lob, (-1)^(s+1))
end

TableauVPLobattoIIIB(s) = TableauVPLobattoIIIB(Float64,s)


@doc raw"""
Tableau for an implicit partitioned Lobatto IIIC Runge-Kutta method with s stages

Uses Lobatto IIIC for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIC(::Type{T}, s) where {T}
    lob = TableauLobattoIIIC(T,s)
    TableauVPRK(Symbol("LobattoIIIC$s"), 2s-2, lob, lob, (-1)^(s+1))
end

TableauVPLobattoIIIC(s) = TableauVPLobattoIIIC(Float64,s)


@doc raw"""
Tableau for an implicit partitioned Lobatto IIIC̄ Runge-Kutta method with s stages

Uses Lobatto IIIC̄ for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIC̄(::Type{T}, s) where {T}
    lob = TableauLobattoIIIC̄(T,s)
    TableauVPRK(Symbol("LobattoIIIC̄$s"), 2s-2, lob, lob, (-1)^(s+1))
end

TableauVPLobattoIIIC̄(s) = TableauVPLobattoIIIC̄(Float64,s)


@doc raw"""
Tableau for variational Lobatto IIID method with s stages

Uses Lobatto IIID for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIID(::Type{T}, s) where {T}
    lob = TableauLobattoIIID(T,s)
    TableauVPRK(Symbol("LobattoIIID$(s)"), 2s-s, lob, lob, (-1)^s)
end

TableauVPLobattoIIID(s) = TableauVPLobattoIIID(Float64,s)


@doc raw"""
Tableau for variational Lobatto IIIE method with s stages

Uses Lobatto IIIE for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIE(::Type{T}, s) where {T}
    lob = TableauLobattoIIIE(T,s)
    TableauVPRK(Symbol("LobattoIIIE$(s)"), 2, lob, lob, (-1)^s)
end

TableauVPLobattoIIIE(s) = TableauVPLobattoIIIE(Float64,s)


@doc raw"""
Tableau for variational Lobatto IIIF method with s stages

Uses Lobatto IIIF for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIF(::Type{T}, s) where {T}
    lob = TableauLobattoIIIF(T,s)
    TableauVPRK(Symbol("LobattoIIIF$(s)"), 2, lob, lob, (-1)^s)
end

TableauVPLobattoIIIF(s) = TableauVPLobattoIIIF(Float64,s)


@doc raw"""
Tableau for variational Lobatto IIIG method with s stages

Uses Lobatto IIIG for both coefficients $a_{ij}$ and $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIG(::Type{T}, s) where {T}
    lob = TableauLobattoIIIG(T,s)
    TableauVPRK(Symbol("LobattoIIIG$(s)"), 2, lob, lob, (-1)^s)
end

TableauVPLobattoIIIG(s) = TableauVPLobattoIIIG(Float64,s)


@doc raw"""
Tableau for variational Lobatto IIIA-IIIB method with s stages

Uses Lobatto IIIA for the coefficients $a_{ij}$ and its symplectic conjugate, Lobatto IIIB, for the coefficients $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIAIIIB(::Type{T}, s) where {T}
    TableauVPRK(Symbol("LobattoIIIA$(s)"), 2s-2, TableauLobattoIIIA(T,s), TableauLobattoIIIB(T,s), (-1)^(s+1), get_lobatto_nullvector(s))
end

TableauVPLobattoIIIAIIIB(s) = TableauVPLobattoIIIAIIIB(Float64,s)


@doc raw"""
Tableau for variational Lobatto IIIB-IIIA method with s stages

Uses Lobatto IIIB for the coefficients $a_{ij}$ and its symplectic conjugate, Lobatto IIIA, for the coefficients $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIBIIIA(::Type{T}, s) where {T}
    TableauVPRK(Symbol("LobattoIIIB$(s)"), 2s-2, TableauLobattoIIIA(T,s), TableauLobattoIIIB(T,s), (-1)^(s+1), get_lobatto_nullvector(s))
end

TableauVPLobattoIIIBIIIA(s) = TableauVPLobattoIIIBIIIA(Float64,s)


@doc raw"""
Tableau for variational Lobatto IIIA method with s stages

Uses Lobatto IIIA for the coefficients $a_{ij}$ and its symplectic conjugate for the coefficients $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIAIIIĀ(::Type{T}, s) where {T}
    TableauVPRK(Symbol("LobattoIIIA$(s)"), 2s-2, TableauLobattoIIIA(T,s), (-1)^(s+1), get_lobatto_nullvector(s))
end

TableauVPLobattoIIIAIIIĀ(s) = TableauVPLobattoIIIAIIIĀ(Float64,s)


@doc raw"""
Tableau for variational Lobatto IIIB method with s stages

Uses Lobatto IIIB for the coefficients $a_{ij}$ and its symplectic conjugate for the coefficients $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIBIIIB̄(::Type{T}, s) where {T}
    TableauVPRK(Symbol("LobattoIIIB$(s)"), 2s-2, TableauLobattoIIIB(T,s), (-1)^(s+1), get_lobatto_nullvector(s))
end

TableauVPLobattoIIIBIIIB̄(s) = TableauVPLobattoIIIBIIIB̄(Float64,s)


@doc raw"""
Tableau for variational Lobatto IIIC method with s stages

Uses Lobatto IIIC for the coefficients $a_{ij}$ and its symplectic conjugate for the coefficients $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIICIIIC̄(::Type{T}, s) where {T}
    TableauVPRK(Symbol("LobattoIIIC$(s)"), 2s-2, TableauLobattoIIIC(T,s), (-1)^(s+1))
end

TableauVPLobattoIIICIIIC̄(s) = TableauVPLobattoIIICIIIC̄(Float64,s)


@doc raw"""
Tableau for variational Lobatto IIIC̄ method with s stages

Uses Lobatto IIIC̄ for the coefficients $a_{ij}$ and its symplectic conjugate for the coefficients $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIC̄IIIC(::Type{T}, s) where {T}
    TableauVPRK(Symbol("LobattoIIIC̄$(s)"), 2s-2, TableauLobattoIIIC̄(T,s), (-1)^(s+1))
end

TableauVPLobattoIIIC̄IIIC(s) = TableauVPLobattoIIIC̄IIIC(Float64,s)


@doc raw"""
Tableau for variational Lobatto IIID method with s stages

Uses Lobatto IIID for the coefficients $a_{ij}$ and its symplectic conjugate for the coefficients $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIDIIID̄(::Type{T}, s) where {T}
    TableauVPRK(Symbol("LobattoIIID$(s)"), 2, TableauLobattoIIIF(T,s), (-1)^s)
end

TableauVPLobattoIIIDIIID̄(s) = TableauVPLobattoIIIDIIID̄(Float64,s)


@doc raw"""
Tableau for variational Lobatto IIIE method with s stages

Uses Lobatto IIIE for the coefficients $a_{ij}$ and its symplectic conjugate for the coefficients $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIEIIIĒ(::Type{T}, s) where {T}
    TableauVPRK(Symbol("LobattoIIIE$(s)"), 2, TableauLobattoIIIF(T,s), (-1)^s)
end

TableauVPLobattoIIIEIIIĒ(s) = TableauVPLobattoIIIEIIIĒ(Float64,s)


@doc raw"""
Tableau for variational Lobatto IIIF method with s stages

Uses Lobatto IIIF for the coefficients $a_{ij}$ and its symplectic conjugate for the coefficients $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIFIIIF̄(::Type{T}, s) where {T}
    TableauVPRK(Symbol("LobattoIIIF$(s)"), 2, TableauLobattoIIIF(T,s), (-1)^s)
end

TableauVPLobattoIIIFIIIF̄(s) = TableauVPLobattoIIIFIIIF̄(Float64,s)


@doc raw"""
Tableau for variational Lobatto IIIF method with s stages

Uses symplectic conjugate of Lobatto IIIF for the coefficients $a_{ij}$ and Lobatto IIIF for the coefficients $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIF̄IIIF(::Type{T}, s) where {T}
    TableauVPRK(Symbol("LobattoIIIF̄$(s)"), 2, TableauLobattoIIIF(T,s), (-1)^s)
end

TableauVPLobattoIIIF̄IIIF(s) = TableauVPLobattoIIIF̄IIIF(Float64,s)


@doc raw"""
Tableau for variational Lobatto IIIG method with s stages

Uses Lobatto IIIG for the coefficients $a_{ij}$ and its symplectic conjugate for the coefficients $\bar{a}_{ij}$.
"""
function TableauVPLobattoIIIGIIIḠ(::Type{T}, s) where {T}
    TableauVPRK(Symbol("LobattoIIIG$(s)"), 2, TableauLobattoIIIF(T,s), (-1)^s)
end

TableauVPLobattoIIIGIIIḠ(s) = TableauVPLobattoIIIGIIIḠ(Float64,s)

