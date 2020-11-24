
function getTableauSLRK(s, o, tsym, q, p, q̃, p̃, ω, d=Nothing)
    TableauSLRK(tsym, o, s, q, p, q̃, p̃, ω, d)
end


function TableauSLRKLobattoIII(s, q, p)
    o = 2s-2
    tsym = Symbol("SLRKLobattoIII")
    getTableauSLRK(s, o, tsym, q, p, q, p, get_lobatto_ω_matrix(s), get_lobatto_d_vector(s))
end

function TableauSLRKLobattoIIIAB(s)
    q = CoefficientsLobattoIIIA(s)
    p = CoefficientsLobattoIIIB(s)
    TableauSLRKLobattoIII(s, q, p)
end

function TableauSLRKLobattoIIIBA(s)
    q = CoefficientsLobattoIIIB(s)
    p = CoefficientsLobattoIIIA(s)
    TableauSLRKLobattoIII(s, q, p)
end

function TableauSLRKLobattoIIICC̄(s)
    q = CoefficientsLobattoIIIC(s)
    p = CoefficientsLobattoIIIC̄(s)
    TableauSLRKLobattoIII(s, q, p)
end

function TableauSLRKLobattoIIIC̄C(s)
    q = CoefficientsLobattoIIIC̄(s)
    p = CoefficientsLobattoIIIC(s)
    TableauSLRKLobattoIII(s, q, p)
end

function TableauSLRKLobattoIIID(s)
    l = CoefficientsLobattoIIID(s)
    TableauSLRKLobattoIII(s, l, l)
end

function TableauSLRKLobattoIIIE(s)
    l = CoefficientsLobattoIIIE(s)
    TableauSLRKLobattoIII(s, l, l)
end
