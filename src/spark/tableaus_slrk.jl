
function SLRKLobattoIII(tsym, s, q, p)
    o = 2s-2
    SLRK(tsym, o, s, q, p, q, p, lobatto_ω_matrix(s), get_lobatto_nullvector(s))
end

function SLRKLobattoIIIAB(s)
    q = TableauLobattoIIIA(s)
    p = TableauLobattoIIIĀ(s)
    SLRKLobattoIII(Symbol("SLRKLobattoIIIAIIIĀ"), s, q, p)
end

function SLRKLobattoIIIBA(s)
    q = TableauLobattoIIIB(s)
    p = TableauLobattoIIIB̄(s)
    SLRKLobattoIII(Symbol("SLRKLobattoIIIBIIIB̄"), s, q, p)
end

function SLRKLobattoIIICC̄(s)
    q = TableauLobattoIIIC(s)
    p = TableauLobattoIIIC̄(s)
    SLRKLobattoIII(Symbol("SLRKLobattoIIICIIIC̄"), s, q, p)
end

function SLRKLobattoIIIC̄C(s)
    q = TableauLobattoIIIC̄(s)
    p = TableauLobattoIIIC(s)
    SLRKLobattoIII(Symbol("SLRKLobattoIIICIIIC̄"), s, q, p)
end

function SLRKLobattoIIID(s)
    q = TableauLobattoIIID(s)
    p = TableauLobattoIIID̄(s)
    SLRKLobattoIII(Symbol("SLRKLobattoIIID"), s, q, p)
end

function SLRKLobattoIIIE(s)
    q = TableauLobattoIIIE(s)
    p = TableauLobattoIIIĒ(s)
    SLRKLobattoIII(Symbol("SLRKLobattoIIIE"), s, q, p)
end
