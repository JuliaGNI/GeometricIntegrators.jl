
function getTableauSLRK(s, o, tsym, q, p, q̃, p̃, ω, d=Nothing)
    TableauSLRK(tsym, o, s, q, p, q̃, p̃, ω, d)
end


function getTableauSLRKLobIII(s, q, p)
    o = 2s-2
    tsym = Symbol("SLRKLobIII")
    getTableauSLRK(s, o, tsym, q, p, q, p, get_lobatto_ω_matrix(s), get_lobatto_d_vector(s))
end

function getTableauSLRKLobIIIAB(s)
    q = getCoefficientsLobIIIA(s)
    p = getCoefficientsLobIIIB(s)
    getTableauSLRKLobIII(s, q, p)
end

function getTableauSLRKLobIIIC(s)
    q = getCoefficientsLobIII(s)
    p = getCoefficientsLobIIIC(s)
    getTableauSLRKLobIII(s, q, p)
end

function getTableauSLRKLobIIID(s)
    l = getCoefficientsLobIIID(s)
    getTableauSLRKLobIII(s, l, l)
end

function getTableauSLRKLobIIIE(s)
    l = getCoefficientsLobIIIE(s)
    getTableauSLRKLobIII(s, l, l)
end
