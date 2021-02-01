
function getTableauSLRK(s, o, tsym, q, p, q̃, p̃, ω, d=Nothing)
    TableauSLRK(tsym, o, s, q, p, q̃, p̃, ω, d)
end


function TableauSLRKLobattoIII(s, q, p)
    o = 2s-2
    tsym = Symbol("SLRKLobattoIII")
    getTableauSLRK(s, o, tsym, q, p, q, p, get_lobatto_ω_matrix(s), get_lobatto_d_vector(s))
end

function TableauSLRKLobattoIIIAB(s)
    q = TableauLobattoIIIA(s)
    p = TableauLobattoIIIB(s)
    TableauSLRKLobattoIII(s, q, p)
end

function TableauSLRKLobattoIIIBA(s)
    q = TableauLobattoIIIB(s)
    p = TableauLobattoIIIA(s)
    TableauSLRKLobattoIII(s, q, p)
end

function TableauSLRKLobattoIIICC̄(s)
    q = TableauLobattoIIIC(s)
    p = TableauLobattoIIIC̄(s)
    TableauSLRKLobattoIII(s, q, p)
end

function TableauSLRKLobattoIIIC̄C(s)
    q = TableauLobattoIIIC̄(s)
    p = TableauLobattoIIIC(s)
    TableauSLRKLobattoIII(s, q, p)
end

function TableauSLRKLobattoIIID(s)
    l = TableauLobattoIIID(s)
    TableauSLRKLobattoIII(s, l, l)
end

function TableauSLRKLobattoIIIE(s)
    l = TableauLobattoIIIE(s)
    TableauSLRKLobattoIII(s, l, l)
end
