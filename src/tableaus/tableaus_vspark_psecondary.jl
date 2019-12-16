
function getTableauVSPARK(s, σ, o, tsym, g, h, lq, lp, ω, d=Nothing)
    # α_q_1 and α_p_2 need to be conjugate symplectic
    # α_q_1 and α_p_3 need to be conjugate symplectic
    # α_q_2 and α_p_2 need to be conjugate symplectic
    # α_q_2 and α_p_3 need to be conjugate symplectic
    # -> α_q_1 and α_q_2 need to be identical
    # -> α_p_2 and α_p_3 need to be identical

    α_q_1 = lq.a
    α_q_2 = lq.a
    α_p_2 = lp.a
    α_p_3 = lp.a

    b_q_1 = lq.b
    b_q_2 = lq.b
    b_p_2 = lp.b
    b_p_3 = lp.b

    # α_p_1 is free, but determines a_q_1 and a_q_2
    α_p_1 = h.a

    # a_p_1, a_p_2 and a_p_3 are free (they are not used in VPSARKpSecondary integrators anymore)
    a_p_1 = g.a
    b_p_1 = g.b

    β_q_1 = lq.b
    β_q_2 = lq.b

    β_p_1 = lp.b
    β_p_2 = lp.b
    β_p_3 = lp.b

    # a_q_1_ij = b_q_1_j / b_p_1_i * (b_p_1_i - α_p_1_ji)
    a_q_1 = zeros(s,σ)
    for i in 1:s
        for j in 1:σ
            a_q_1[i,j] = b_q_1[j] / b_p_1[i] * (b_p_1[i] - α_p_1[j,i])
        end
    end

    # a_q_2_ij = b_q_2_j / b_p_1_i * (b_p_1_i - α_p_1_ji)
    a_q_2 = zeros(s,σ)
    for i in 1:s
        for j in 1:σ
            a_q_2[i,j] = b_q_2[j] / b_p_1[i] * (b_p_1[i] - α_p_1[j,i])
        end
    end

    a_p_2 = a_q_1
    a_p_3 = a_q_2


    a_q = (a_q_1, a_q_2)
    b_q = (b_q_1, b_q_2)
    c_q = g.c

    a_p = (a_p_1, a_p_2, a_p_3)
    b_p = (b_p_1, b_p_2, b_p_3)
    c_p = g.c

    α_q = (α_q_1, α_q_2)
    β_q = (β_q_1, β_q_2)
    γ_q = lq.c

    α_p = (α_p_1, α_p_2, α_p_3)
    β_p = (β_p_1, β_p_2, β_p_3)
    γ_p = lp.c

    coeff_q = CoefficientsSPARK(tsym, o, s, σ, a_q, b_q, c_q)
    coeff_p = CoefficientsSPARK(tsym, o, s, σ, a_p, b_p, c_p)
    coeff_q̃ = CoefficientsSPARK(tsym, o, σ, s, α_q, β_q, γ_q)
    coeff_p̃ = CoefficientsSPARK(tsym, o, σ, s, α_p, β_p, γ_p)

    TableauVSPARKsecondary(tsym, o, s, σ, coeff_q, coeff_p, coeff_q̃, coeff_p̃, ω, d)
end


function getTableauVSPARKLobIII(s, lq, lp)
    o = 2s-2
    tsym = Symbol("VSPARKLobIII")
    getTableauVSPARK(s, s, o, tsym, lq, lp, lq, lp, get_lobatto_ω_matrix(s), get_lobatto_d_vector(s))
end

function getTableauVSPARKLobIIIAB(s)
    lq = getCoefficientsLobIIIA(s)
    lp = getCoefficientsLobIIIB(s)
    getTableauVSPARKLobIII(s, lq, lp)
end

function getTableauVSPARKLobIIIC(s)
    l = getCoefficientsLobIII(s)
    getTableauVSPARKLobIII(s, l, l)
end

function getTableauVSPARKLobIIID(s)
    l = getCoefficientsLobIIID(s)
    getTableauVSPARKLobIII(s, l, l)
end

function getTableauVSPARKLobIIIE(s)
    l = getCoefficientsLobIIIE(s)
    getTableauVSPARKLobIII(s, l, l)
end


function getTableauVSPARKGLRKLobIII(s, σ, lq, lp)
    o = 2s
    g = getCoefficientsGLRK(s)
    tsym = Symbol("VSPARKGLRKLobIII")
    getTableauVSPARK(s, σ, o, tsym, g, get_lobatto_interstage_coefficients(s, σ), lq, lp, get_GLRK_ω_matrix(σ), get_lobatto_d_vector(σ))
end

function getTableauVSPARKGLRKLobIIIAB(s, σ=s+1)
    getTableauVSPARKGLRKLobIII(s, σ, getCoefficientsLobIIIA(σ), getCoefficientsLobIIIB(σ))
end

function getTableauVSPARKGLRKLobIIIC(s, σ=s+1)
    getTableauVSPARKGLRKLobIII(s, σ, getCoefficientsLobIII(σ), getCoefficientsLobIIIC(σ))
end

function getTableauVSPARKGLRKLobIIID(s, σ=s+1)
    getTableauVSPARKGLRKLobIII(s, σ, getCoefficientsLobIIID(σ), getCoefficientsLobIIID(σ))
end

function getTableauVSPARKGLRKLobIIIE(s, σ=s+1)
    getTableauVSPARKGLRKLobIII(s, σ, getCoefficientsLobIIIE(σ), getCoefficientsLobIIIE(σ))
end
