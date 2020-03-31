
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


function getTableauHSPARKLobIII(s, lq, lp)
    o = 2s-2
    tsym = Symbol("HSPARKLobIII")
    getTableauHSPARK(s, s, o, tsym, lq, lp, lq, lp, get_lobatto_ω_matrix(s), get_lobatto_d_vector(s))
end

function getTableauHSPARKLobIIIAB(s)
    lq = getCoefficientsLobIIIA(s)
    lp = getCoefficientsLobIIIB(s)
    getTableauHSPARKLobIII(s, lq, lp)
end

function getTableauHSPARKLobIIIC(s)
    lq = getCoefficientsLobIII(s)
    lp = getCoefficientsLobIIIC(s)
    getTableauHSPARKLobIII(s, lq, lp)
end

function getTableauHSPARKLobIIID(s)
    l = getCoefficientsLobIIID(s)
    getTableauHSPARKLobIII(s, l, l)
end

function getTableauHSPARKLobIIIE(s)
    l = getCoefficientsLobIIIE(s)
    getTableauHSPARKLobIII(s, l, l)
end


function getTableauHSPARKGLRKLobIII(s, σ, lq, lp)
    o = 2s
    g = getCoefficientsGLRK(s)
    tsym = Symbol("HSPARKGLRKLobIII")
    getTableauHSPARK(s, σ, o, tsym, g, get_lobatto_interstage_coefficients(s, σ), lq, lp, get_GLRK_ω_matrix(σ), get_lobatto_d_vector(σ))
end

function getTableauHSPARKGLRKLobIIIAB(s, σ=s+1)
    getTableauHSPARKGLRKLobIII(s, σ, getCoefficientsLobIIIA(σ), getCoefficientsLobIIIB(σ))
end

function getTableauHSPARKGLRKLobIIIC(s, σ=s+1)
    getTableauHSPARKGLRKLobIII(s, σ, getCoefficientsLobIII(σ), getCoefficientsLobIIIC(σ))
end

function getTableauHSPARKGLRKLobIIID(s, σ=s+1)
    getTableauHSPARKGLRKLobIII(s, σ, getCoefficientsLobIIID(σ), getCoefficientsLobIIID(σ))
end

function getTableauHSPARKGLRKLobIIIE(s, σ=s+1)
    getTableauHSPARKGLRKLobIII(s, σ, getCoefficientsLobIIIE(σ), getCoefficientsLobIIIE(σ))
end
