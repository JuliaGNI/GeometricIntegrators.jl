
function getTableauHSPARK(s, σ, o, tsym, g, h, lq, lp, ω, d=Nothing)
    # α_q_1 and α_p_2 need to be conjugate symplectic
    # α_q_1 and α_p_3 need to be conjugate symplectic
    # α_q_2 and α_p_2 need to be conjugate symplectic
    # α_q_2 and α_p_3 need to be conjugate symplectic
    # -> α_q_1 and α_q_2 need to be identical
    # -> α_p_2 and α_p_3 need to be identical

    α_q_2 = lq.a
    α_q_3 = lq.a
    α_p_2 = lp.a
    α_p_3 = lp.a

    b_q_2 = lq.b
    b_q_3 = lq.b
    b_p_2 = lp.b
    b_p_3 = lp.b

    # α_p_1 is free, but determines a_q_1 and a_q_2
    α_q_1 = h.a
    α_p_1 = h.a

    # a_p_1, a_p_2 and a_p_3 are free (they are not used in HSPARKpSecondary integrators anymore)
    # REALLY?
    a_q_1 = g.a
    b_q_1 = g.b

    a_p_1 = g.a
    b_p_1 = g.b

    a_p_2 = g.a
    a_p_3 = g.a

    β_q_1 = lq.b
    β_q_2 = lq.b
    β_q_3 = lq.b

    β_p_1 = lp.b
    β_p_2 = lp.b
    β_p_3 = lp.b

    # a_q_1_ij = b_q_1_j / b_p_1_i * (b_p_1_i - α_p_1_ji)
    # a_q_1 = zeros(s,σ)
    # for i in 1:s
    #     for j in 1:σ
    #         a_q_1[i,j] = b_q_1[j] / b_p_1[i] * (b_p_1[i] - α_p_1[j,i])
    #     end
    # end

    # a_q_1_ij = b_q_1_j / b_p_1_i * (b_p_1_i - α_p_1_ji)
    # TODO check this
    a_q_2 = zeros(s,σ)
    for i in 1:s
        for j in 1:σ
            a_q_2[i,j] = b_q_2[j] / b_p_1[i] * (b_p_1[i] - α_p_2[j,i])
        end
    end

    # a_q_2_ij = b_q_2_j / b_p_1_i * (b_p_1_i - α_p_1_ji)
    # TODO check this
    a_q_3 = zeros(s,σ)
    for i in 1:s
        for j in 1:σ
            a_q_3[i,j] = b_q_3[j] / b_p_1[i] * (b_p_1[i] - α_p_3[j,i])
        end
    end


    a_q = (a_q_1, a_q_2, a_q_3)
    b_q = (b_q_1, b_q_2, b_q_3)
    c_q = g.c

    a_p = (a_p_1, a_p_2, a_p_3)
    b_p = (b_p_1, b_p_2, b_p_3)
    c_p = g.c

    α_q = (α_q_1, α_q_2, α_q_3)
    β_q = (β_q_1, β_q_2, β_q_3)
    γ_q = lq.c

    α_p = (α_p_1, α_p_2, α_p_3)
    β_p = (β_p_1, β_p_2, β_p_3)
    γ_p = lp.c

    coeff_q = CoefficientsSPARK(tsym, o, s, σ, a_q, b_q, c_q)
    coeff_p = CoefficientsSPARK(tsym, o, s, σ, a_p, b_p, c_p)
    coeff_q̃ = CoefficientsSPARK(tsym, o, σ, s, α_q, β_q, γ_q)
    coeff_p̃ = CoefficientsSPARK(tsym, o, σ, s, α_p, β_p, γ_p)

    TableauHSPARKsecondary(tsym, o, s, σ, coeff_q, coeff_p, coeff_q̃, coeff_p̃, ω, d)
end


function TableauHSPARKLobattoIII(s, lq, lp; name = Symbol("HSPARKLobattoIII"))
    o = 2s-2
    getTableauHSPARK(s, s, o, name, lq, lp, lq, lp, get_lobatto_ω_matrix(s), get_lobatto_nullvector(s))
end

function TableauHSPARKLobattoIIIAB(s)
    lq = TableauLobattoIIIA(s)
    lp = TableauLobattoIIIB(s)
    TableauHSPARKLobattoIII(s, lq, lp; name = Symbol("HSPARKLobattoIIIAIIIB"))
end

function TableauHSPARKLobattoIIIBA(s)
    lq = TableauLobattoIIIB(s)
    lp = TableauLobattoIIIA(s)
    TableauHSPARKLobattoIII(s, lq, lp; name = Symbol("HSPARKLobattoIIIBIIIA"))
end

function TableauHSPARKLobattoIIICC̄(s)
    lq = TableauLobattoIIIC(s)
    lp = TableauLobattoIIIC̄(s)
    TableauHSPARKLobattoIII(s, lq, lp; name = Symbol("HSPARKLobattoIIICIIIC̄"))
end

function TableauHSPARKLobattoIIIC̄C(s)
    lq = TableauLobattoIIIC̄(s)
    lp = TableauLobattoIIIC(s)
    TableauHSPARKLobattoIII(s, lq, lp; name = Symbol("HSPARKLobattoIIIC̄IIIC"))
end

function TableauHSPARKLobattoIIID(s)
    l = TableauLobattoIIID(s)
    TableauHSPARKLobattoIII(s, l, l; name = Symbol("HSPARKLobattoIIID"))
end

function TableauHSPARKLobattoIIIE(s)
    l = TableauLobattoIIIE(s)
    TableauHSPARKLobattoIII(s, l, l; name = Symbol("HSPARKLobattoIIIE"))
end


function TableauHSPARKGLRKLobattoIII(s, σ, lq, lp; name = Symbol("HSPARKGLRKLobattoIII"))
    o = 2s
    g = TableauGauss(s)
    getTableauHSPARK(s, σ, o, name, g, get_lobatto_glrk_coefficients(s, σ), lq, lp, get_GLRK_ω_matrix(σ), get_lobatto_nullvector(σ))
end

function TableauHSPARKGLRKLobattoIIIAB(s, σ=s+1)
    TableauHSPARKGLRKLobattoIII(s, σ, TableauLobattoIIIA(σ), TableauLobattoIIIB(σ); name = Symbol("HSPARKGLRKLobattoIIIAIIIB"))
end

function TableauHSPARKGLRKLobattoIIIBA(s, σ=s+1)
    TableauHSPARKGLRKLobattoIII(s, σ, TableauLobattoIIIB(σ), TableauLobattoIIIA(σ); name = Symbol("HSPARKGLRKLobattoIIIBIIIA"))
end

function TableauHSPARKGLRKLobattoIIICC̄(s, σ=s+1)
    TableauHSPARKGLRKLobattoIII(s, σ, TableauLobattoIIIC(σ), TableauLobattoIIIC̄(σ); name = Symbol("HSPARKGLRKLobattoIIICIIIC̄"))
end

function TableauHSPARKGLRKLobattoIIIC̄C(s, σ=s+1)
    TableauHSPARKGLRKLobattoIII(s, σ, TableauLobattoIIIC̄(σ), TableauLobattoIIIC(σ); name = Symbol("HSPARKGLRKLobattoIIIC̄IIIC"))
end

function TableauHSPARKGLRKLobattoIIID(s, σ=s+1)
    TableauHSPARKGLRKLobattoIII(s, σ, TableauLobattoIIID(σ), TableauLobattoIIID(σ); name = Symbol("HSPARKGLRKLobattoIIID"))
end

function TableauHSPARKGLRKLobattoIIIE(s, σ=s+1)
    TableauHSPARKGLRKLobattoIII(s, σ, TableauLobattoIIIE(σ), TableauLobattoIIIE(σ); name = Symbol("HSPARKGLRKLobattoIIIE"))
end
