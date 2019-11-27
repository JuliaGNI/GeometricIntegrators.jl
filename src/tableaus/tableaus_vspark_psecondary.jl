
function getTableauVSPARKGLRK(s, σ=s)
    o = 2s

    glrk_q = getCoefficientsGLRK(s)
    glrk_p = get_symplectic_conjugate_coefficients(glrk_q)

    # a_q_1 and a_p_2 need to be conjugate symplectic
    a_q_1 = glrk_q.a
    b_q_1 = glrk_q.b

    a_p_1 = glrk_p.a
    b_p_1 = glrk_p.b

    a_p_2 = glrk_p.a
    b_p_2 = glrk_p.b

    # α_q_2 and α_p_3 need to be conjugate symplectic
    # α_q_1 and α_p_2 are free

    α_q_2 = glrk_q.a
    α_p_3 = glrk_p.a
    α_q_1 = glrk_q.a
    α_p_2 = glrk_p.a

    b_q_2 = glrk_q.b
    b_p_3 = glrk_p.b

    β_q_1 = glrk_q.b
    β_q_2 = glrk_q.b

    β_p_1 = glrk_p.b
    β_p_2 = glrk_p.b
    β_p_3 = glrk_p.b

    # a_p_3_ij = b_p_3_j / b_q_1_i * (b_q_1_i - α_q_1_ji)
    a_p_3 = zeros(s,σ)
    for i in 1:s
        for j in 1:σ
            a_p_3[i,j] = b_p_3[j] / b_q_1[i] * (b_q_1[i] - α_q_1[j,i])
        end
    end

    # a_q_2_ij = b_q_2_j / b_p_2_i * (b_p_2_i - α_p_2_ji)
    a_q_2 = zeros(s,σ)
    for i in 1:s
        for j in 1:σ
            a_q_2[i,j] = b_q_2[j] / b_p_2[i] * (b_p_2[i] - α_p_2[j,i])
        end
    end

    # α_p_1_ij = b_p_1_j / b_q_2_i * (b_q_2_i - a_q_2_ji)
    α_p_1 = zeros(σ,s)
    for i in 1:σ
        for j in 1:s
            α_p_1[i,j] = b_p_1[j] / b_q_2[i] * (b_q_2[i] - a_q_2[j,i])
        end
    end

    a_q = (a_q_1, a_q_2)
    b_q = (b_q_1, b_q_2)
    c_q = glrk_q.c

    a_p = (a_p_1, a_p_2, a_p_3)
    b_p = (b_p_1, b_p_2, b_p_3)
    c_p = glrk_p.c

    α_q = (α_q_1, α_q_2)
    β_q = (β_q_1, β_q_2)
    γ_q = glrk_q.c

    α_p = (α_p_1, α_p_2, α_p_3)
    β_p = (β_p_1, β_p_2, β_p_3)
    γ_p = glrk_p.c

    coeff_q = CoefficientsSPARK(Symbol("VSPARKGLRK"),   o, s, σ, a_q, b_q, c_q)
    coeff_p = CoefficientsSPARK(Symbol("VSPARKGLRK"),   o, s, σ, a_p, b_p, c_p)
    coeff_q̃ = CoefficientsSPARK(Symbol("VSPARKGLRK"), o, σ, s, α_q, β_q, γ_q)
    coeff_p̃ = CoefficientsSPARK(Symbol("VSPARKGLRK"), o, σ, s, α_p, β_p, γ_p)

    ω = zeros(σ, σ+1)
    g = glrk_q
    ω[σ,σ+1] = 1
    for i in 1:σ-1
        for j in 1:σ
            ω[i,j] = g.b[j] * g.c[j]^(i-1)
        end
    end

    TableauVSPARKsecondary(Symbol("VSPARKGLRK"), o, s, σ, coeff_q, coeff_p, coeff_q̃, coeff_p̃, ω)
end



function get_lobatto_d_vector(σ)
    if σ == 2
        d = [+1.0, -1.0]
    elseif σ == 3
        d = [+1.0, -2.0, +1.0]
    elseif σ == 4
        d = [+1.0, -√5, +√5, -1.0]
    else
        @error("We don't have a d vector for σ=" * str(σ) * ".")
    end
    return d
end


function getTableauVSPARKLobIII(s)
    σ = s
    o = 2s-2

    lobIII  = getCoefficientsLobIII(σ)
    lobIIIA = getCoefficientsLobIIIA(σ)
    lobIIIB = getCoefficientsLobIIIB(σ)
    lobIIIC = getCoefficientsLobIIIC(σ)
    lobIIID = getCoefficientsLobIIID(σ)
    lobIIIE = getCoefficientsLobIIIE(σ)

    # α_q_1 and α_p_2 need to be conjugate symplectic
    # α_q_1 and α_p_3 need to be conjugate symplectic
    # α_q_2 and α_p_2 need to be conjugate symplectic
    # α_q_2 and α_p_3 need to be conjugate symplectic
    # -> α_q_1 and α_q_2 need to be identical
    # -> α_p_2 and α_p_3 need to be identical

    α_q_1 = lobIIIE.a
    α_q_2 = lobIIIE.a
    α_p_2 = lobIIIE.a
    α_p_3 = lobIIIE.a

    b_q_1 = lobIIIE.b
    b_q_2 = lobIIIE.b
    b_p_2 = lobIIIE.b
    b_p_3 = lobIIIE.b

    # α_p_1 is free, but determines a_q_1 and a_q_2

    α_p_1 = lobIIIB.a#[1:σ,1:σ-1].a
    b_p_1 = lobIIIB.b

    # a_p_1, a_p_2 and a_p_3 are free

    a_p_1 = lobIIID.a
    b_p_1 = lobIIID.b

    # a_p_2 = lobIIIA.a[2:σ,1:σ]
    # a_p_3 = lobIIIA.a[2:σ,1:σ]
    a_p_2 = lobIIID.a
    a_p_3 = lobIIID.a


    β_q_1 = lobIIIA.b
    β_q_2 = lobIIIA.b

    β_p_1 = lobIIID.b
    β_p_2 = lobIIIB.b
    β_p_3 = lobIIIB.b


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


    a_q = (a_q_1, a_q_2)
    b_q = (b_q_1, b_q_2)
    c_q = lobIII.c

    a_p = (a_p_1, a_p_2, a_p_3)
    b_p = (b_p_1, b_p_2, b_p_3)
    c_p = lobIII.c

    α_q = (α_q_1, α_q_2)
    β_q = (β_q_1, β_q_2)
    γ_q = lobIII.c

    α_p = (α_p_1, α_p_2, α_p_3)
    β_p = (β_p_1, β_p_2, β_p_3)
    γ_p = lobIII.c

    coeff_q = CoefficientsSPARK(Symbol("VSPARKLobIII"),   o, s, σ, a_q, b_q, c_q)
    coeff_p = CoefficientsSPARK(Symbol("VSPARKLobIII"),   o, s, σ, a_p, b_p, c_p)
    coeff_q̃ = CoefficientsSPARK(Symbol("VSPARKLobIII"), o, σ, s, α_q, β_q, γ_q)
    coeff_p̃ = CoefficientsSPARK(Symbol("VSPARKLobIII"), o, σ, s, α_p, β_p, γ_p)

    ω = zeros(σ, σ+1)
    ω[1:σ-1,1:σ] .= lobIIIA.a[2:σ,1:σ]
    ω[σ,σ+1] = 1

    TableauVSPARKsecondary(Symbol("VSPARKLobIII"), o, s, σ, coeff_q, coeff_p, coeff_q̃, coeff_p̃, ω)#, get_lobatto_d_vector(σ))
end


function getTableauVSPARKGLRKLobIII(s, σ=s+1)
    o = 2s

    glrk    = getCoefficientsGLRK(s)
    lobIII  = getCoefficientsLobIII(σ)
    lobIIIA = getCoefficientsLobIIIA(σ)
    lobIIIB = getCoefficientsLobIIIB(σ)
    lobIIIC = getCoefficientsLobIIIC(σ)
    lobIIID = getCoefficientsLobIIID(σ)
    lobIIIE = getCoefficientsLobIIIE(σ)

    # α_q_1 and α_p_2 need to be conjugate symplectic
    # α_q_1 and α_p_3 need to be conjugate symplectic
    # α_q_2 and α_p_2 need to be conjugate symplectic
    # α_q_2 and α_p_3 need to be conjugate symplectic
    # -> α_q_1 and α_q_2 need to be identical
    # -> α_p_2 and α_p_3 need to be identical

    α_q_1 = lobIIIA.a
    α_q_2 = lobIIIA.a
    α_p_2 = lobIIIB.a
    α_p_3 = lobIIIB.a

    b_q_1 = lobIIIA.b
    b_q_2 = lobIIIA.b
    b_p_2 = lobIIIB.b
    b_p_3 = lobIIIB.b

    # α_p_1 is free, but determines a_q_1 and a_q_2

    if σ == 2
        α_p_1  = zeros(Float64, 2, 1)
        α_p_1 .= [[0//1]
                  [1//1]]
    elseif σ == 3
        α_p_1  = zeros(Float64, 3, 2)
        α_p_1 .= [[0            0          ]
                  [1//4 + √3/8  1//4 - √3/8]
                  [1//2         1//2       ]]
    else
        @error("σ is not supported")
    end

    # a_p_1, a_p_2 and a_p_3 are free

    a_p_1 = glrk.a
    b_p_1 = glrk.b


    β_q_1 = lobIII.b
    β_q_2 = lobIII.b

    β_p_1 = lobIII.b
    β_p_2 = lobIII.b
    β_p_3 = lobIII.b


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
    c_q = glrk.c

    a_p = (a_p_1, a_p_2, a_p_3)
    b_p = (b_p_1, b_p_2, b_p_3)
    c_p = glrk.c

    α_q = (α_q_1, α_q_2)
    β_q = (β_q_1, β_q_2)
    γ_q = lobIII.c

    α_p = (α_p_1, α_p_2, α_p_3)
    β_p = (β_p_1, β_p_2, β_p_3)
    γ_p = lobIII.c

    coeff_q = CoefficientsSPARK(Symbol("VSPARKGLRKLobIII"), o, s, σ, a_q, b_q, c_q)
    coeff_p = CoefficientsSPARK(Symbol("VSPARKGLRKLobIII"), o, s, σ, a_p, b_p, c_p)
    coeff_q̃ = CoefficientsSPARK(Symbol("VSPARKGLRKLobIII"), o, σ, s, α_q, β_q, γ_q)
    coeff_p̃ = CoefficientsSPARK(Symbol("VSPARKGLRKLobIII"), o, σ, s, α_p, β_p, γ_p)

    ω = zeros(σ, σ+1)
    ω[1:σ-1,1:σ] .= lobIIIA.a[2:σ,1:σ]
    ω[σ,σ+1] = 1

    TableauVSPARKsecondary(Symbol("VSPARKGLRKLobIII"), o, s, σ, coeff_q, coeff_p, coeff_q̃, coeff_p̃, ω)#, get_lobatto_d_vector(σ))
end


# function getTableauVSPARKGLRKLobIII(s, σ=s+1)
#     o = 2s
#
#     glrk_q = getCoefficientsGLRK(s)
#     glrk_p = get_symplectic_conjugate_coefficients(glrk_q)
#
#     # a_q_1 and a_p_2 need to be conjugate symplectic
#     a_q_1 = glrk_q.a
#     b_q_1 = glrk_q.b
#
#     a_p_1 = glrk_p.a
#     b_p_1 = glrk_p.b
#
#     a_p_2 = glrk_p.a
#     b_p_2 = glrk_p.b
#
#     lobIII  = getCoefficientsLobIII(σ)
#     lobIIIA = getCoefficientsLobIIIA(σ)
#     lobIIIB = getCoefficientsLobIIIB(σ)
#     lobIIIC = getCoefficientsLobIIIC(σ)
#     lobIIID = getCoefficientsLobIIID(σ)
#     lobIIIE = getCoefficientsLobIIIE(σ)
#
#     # α_q_2 and α_p_3 need to be conjugate symplectic
#     # α_q_1 and α_p_2 are free
#
#     α_q_2 = lobIIIA.a
#     α_p_3 = lobIIIB.a
#     α_q_1 = lobIIID.a
#     α_p_2 = lobIIID.a
#
#     # α_q_2 = lobIIID.a
#     # α_p_3 = lobIIID.a
#     # α_q_1 = lobIIIE.a
#     # α_p_2 = lobIIIE.a
#
#     # α_q_2 = lobIIID.a
#     # α_p_3 = lobIIID.a
#     # α_q_1 = lobIIIE.a
#     # α_p_2 = lobIIIE.a
#
#     b_q_2 = lobIIIA.b
#     b_p_3 = lobIIIB.b
#
#     β_q_1 = glrk_q.b
#     β_q_2 = lobIIIA.b
#
#     β_p_1 = glrk_p.b
#     β_p_2 = lobIIID.b
#     β_p_3 = lobIIIE.b
#
#     # a_p_3_ij = b_p_3_j / b_q_1_i * (b_q_1_i - α_q_1_ji)
#     a_p_3 = zeros(s,σ)
#     for i in 1:s
#         for j in 1:σ
#             a_p_3[i,j] = b_p_3[j] / b_q_1[i] * (b_q_1[i] - α_q_1[j,i])
#         end
#     end
#
#     # a_q_2_ij = b_q_2_j / b_p_2_i * (b_p_2_i - α_p_2_ji)
#     a_q_2 = zeros(s,σ)
#     for i in 1:s
#         for j in 1:σ
#             a_q_2[i,j] = b_q_2[j] / b_p_2[i] * (b_p_2[i] - α_p_2[j,i])
#         end
#     end
#
#     # α_p_1_ij = b_p_1_j / b_q_2_i * (b_q_2_i - a_q_2_ji)
#     α_p_1 = zeros(σ,s)
#     for i in 1:σ
#         for j in 1:s
#             α_p_1[i,j] = b_p_1[j] / b_q_2[i] * (b_q_2[i] - a_q_2[j,i])
#         end
#     end
#
#     a_q = (a_q_1, a_q_2)
#     b_q = (b_q_1, b_q_2)
#     c_q = glrk_q.c
#
#     a_p = (a_p_1, a_p_2, a_p_3)
#     b_p = (b_p_1, b_p_2, b_p_3)
#     c_p = glrk_p.c
#
#     α_q = (α_q_1, α_q_2)
#     β_q = (β_q_1, β_q_2)
#     γ_q = lobIII.c
#
#     α_p = (α_p_1, α_p_2, α_p_3)
#     β_p = (β_p_1, β_p_2, β_p_3)
#     γ_p = lobIII.c
#
#     coeff_q = CoefficientsSPARK(Symbol("VSPARKGLRK"),   o, s, σ, a_q, b_q, c_q)
#     coeff_p = CoefficientsSPARK(Symbol("VSPARKGLRK"),   o, s, σ, a_p, b_p, c_p)
#     coeff_q̃ = CoefficientsSPARK(Symbol("VSPARKLobIII"), o, σ, s, α_q, β_q, γ_q)
#     coeff_p̃ = CoefficientsSPARK(Symbol("VSPARKLobIII"), o, σ, s, α_p, β_p, γ_p)
#
#     # ω = zeros(σ, σ+1)
#     # ω[1:σ-1,1:σ] .= lobIIIA.a[2:σ,1:σ]
#     # ω[σ,σ+1] = 1
#
#     ω = zeros(σ, σ+1)
#     g = getCoefficientsGLRK(σ)
#     ω[σ,σ+1] = 1
#     for i in 1:σ-1
#         for j in 1:σ
#             ω[i,j] = g.b[j] * g.c[j]^(i-1)
#         end
#     end
#
#     TableauVSPARKsecondary(Symbol("VSPARKGLRKLobIII"), o, s, σ, coeff_q, coeff_p, coeff_q̃, coeff_p̃, ω)#, get_lobatto_d_vector(σ))
# end


function getTableauVSPARKGLRKLobIIIAB(s, σ=s+1)
    o = 2s

    glrk_q = getCoefficientsGLRK(s)
    glrk_p = get_symplectic_conjugate_coefficients(glrk_q)

    # a_q_1 and a_p_2 need to be conjugate symplectic
    a_q_1 = glrk_q.a
    b_q_1 = glrk_q.b

    a_p_1 = glrk_p.a
    b_p_1 = glrk_p.b

    a_p_2 = glrk_p.a
    b_p_2 = glrk_p.b

    lobIIIA = getCoefficientsLobIIIA(σ)
    lobIIIB = getCoefficientsLobIIIB(σ)

    # α_q_2 and α_p_3 need to be conjugate symplectic
    # α_q_1 and α_p_2 are free

    α_q_2 = lobIIIA.a
    α_p_3 = lobIIIB.a
    α_q_1 = lobIIIA.a
    α_p_2 = lobIIIB.a

    b_q_2 = lobIIIA.b
    b_p_3 = lobIIIB.b

    β_q_1 = glrk_q.b
    β_q_2 = lobIIIA.b

    β_p_1 = glrk_p.b
    β_p_2 = lobIIIB.b
    β_p_3 = lobIIIB.b

    # a_p_3_ij = b_p_3_j / b_q_1_i * (b_q_1_i - α_q_1_ji)
    a_p_3 = zeros(s,σ)
    for i in 1:s
        for j in 1:σ
            a_p_3[i,j] = b_p_3[j] / b_q_1[i] * (b_q_1[i] - α_q_1[j,i])
        end
    end

    # a_q_2_ij = b_q_2_j / b_p_2_i * (b_p_2_i - α_p_2_ji)
    a_q_2 = zeros(s,σ)
    for i in 1:s
        for j in 1:σ
            a_q_2[i,j] = b_q_2[j] / b_p_2[i] * (b_p_2[i] - α_p_2[j,i])
        end
    end

    # α_p_1_ij = b_p_1_j / b_q_2_i * (b_q_2_i - a_q_2_ji)
    α_p_1 = zeros(σ,s)
    for i in 1:σ
        for j in 1:s
            α_p_1[i,j] = b_p_1[j] / b_q_2[i] * (b_q_2[i] - a_q_2[j,i])
        end
    end

    a_q = (a_q_1, a_q_2)
    b_q = (b_q_1, b_q_2)
    c_q = glrk_q.c

    a_p = (a_p_1, a_p_2, a_p_3)
    b_p = (b_p_1, b_p_2, b_p_3)
    c_p = glrk_p.c

    α_q = (α_q_1, α_q_2)
    β_q = (β_q_1, β_q_2)
    γ_q = lobIIIB.c

    α_p = (α_p_1, α_p_2, α_p_3)
    β_p = (β_p_1, β_p_2, β_p_3)
    γ_p = lobIIIB.c

    coeff_q = CoefficientsSPARK(Symbol("VSPARKGLRK"),   o, s, σ, a_q, b_q, c_q)
    coeff_p = CoefficientsSPARK(Symbol("VSPARKGLRK"),   o, s, σ, a_p, b_p, c_p)
    coeff_q̃ = CoefficientsSPARK(Symbol("VSPARKLobIII"), o, σ, s, α_q, β_q, γ_q)
    coeff_p̃ = CoefficientsSPARK(Symbol("VSPARKLobIII"), o, σ, s, α_p, β_p, γ_p)

    ω = zeros(σ, σ+1)
    ω[1:σ-1,1:σ] .= lobIIIA.a[2:σ,1:σ]
    ω[σ,σ+1] = 1

    # ω = zeros(σ, σ+1)
    # g = getCoefficientsGLRK(σ)
    # ω[σ,σ+1] = 1
    # for i in 1:σ-1
    #     for j in 1:σ
    #         ω[i,j] = g.b[j] * g.c[j]^(i-1)
    #     end
    # end

    TableauVSPARKsecondary(Symbol("VSPARKGLRKLobIIIAB"), o, s, σ, coeff_q, coeff_p, coeff_q̃, coeff_p̃, ω, get_lobatto_d_vector(σ))
end


function getTableauVSPARKGLRKLobIIIAC(s, σ=s+1)
    o = 2s

    glrk_q = getCoefficientsGLRK(s)
    glrk_p = get_symplectic_conjugate_coefficients(glrk_q)

    a_q_1 = glrk_q.a
    b_q_1 = glrk_q.b

    a_p_1 = glrk_p.a
    b_p_1 = glrk_p.b

    a_p_2 = glrk_p.a
    b_p_2 = glrk_p.b

    lobIII  = getCoefficientsLobIII(σ)
    lobIIIA = getCoefficientsLobIIIA(σ)
    lobIIIC = getCoefficientsLobIIIC(σ)

    # α_q_2 and α_p_3 need to be conjugate symplectic
    # α_q_1 and α_p_2 are free

    α_q_2 = lobIIIC.a
    α_p_3 = lobIII.a
    α_q_1 = lobIIIC.a
    α_p_2 = lobIII.a

    b_q_2 = lobIIIC.b
    b_p_3 = lobIII.b

    β_q_1 = glrk_q.b
    β_q_2 = lobIII.b

    β_p_1 = glrk_p.b
    β_p_2 = lobIII.b
    β_p_3 = lobIII.b

    # a_p_3_ij = b_p_3_j / b_q_1_i * (b_q_1_i - α_q_1_ji)
    a_p_3 = zeros(s,σ)
    for i in 1:s
        for j in 1:σ
            a_p_3[i,j] = b_p_3[j] / b_q_1[i] * (b_q_1[i] - α_q_1[j,i])
        end
    end

    # a_q_2_ij = b_q_2_j / b_p_2_i * (b_p_2_i - α_p_2_ji)
    a_q_2 = zeros(s,σ)
    for i in 1:s
        for j in 1:σ
            a_q_2[i,j] = b_q_2[j] / b_p_2[i] * (b_p_2[i] - α_p_2[j,i])
        end
    end

    # α_p_1_ij = b_p_1_j / b_q_2_i * (b_q_2_i - a_q_2_ji)
    α_p_1 = zeros(σ,s)
    for i in 1:σ
        for j in 1:s
            α_p_1[i,j] = b_p_1[j] / b_q_2[i] * (b_q_2[i] - a_q_2[j,i])
        end
    end

    a_q = (a_q_1, a_q_2)
    b_q = (b_q_1, b_q_2)
    c_q = glrk_q.c

    a_p = (a_p_1, a_p_2, a_p_3)
    b_p = (b_p_1, b_p_2, b_p_3)
    c_p = glrk_p.c

    α_q = (α_q_1, α_q_2)
    β_q = (β_q_1, β_q_2)
    γ_q = lobIIIC.c

    α_p = (α_p_1, α_p_2, α_p_3)
    β_p = (β_p_1, β_p_2, β_p_3)
    γ_p = lobIIIC.c

    coeff_q = CoefficientsSPARK(Symbol("VSPARKGLRK"),   o, s, σ, a_q, b_q, c_q)
    coeff_p = CoefficientsSPARK(Symbol("VSPARKGLRK"),   o, s, σ, a_p, b_p, c_p)
    coeff_q̃ = CoefficientsSPARK(Symbol("VSPARKLobIII"), o, σ, s, α_q, β_q, γ_q)
    coeff_p̃ = CoefficientsSPARK(Symbol("VSPARKLobIII"), o, σ, s, α_p, β_p, γ_p)

    # ω = zeros(σ, σ+1)
    # ω[1:σ-1,1:σ] .= lobIIIA.a[2:σ,1:σ]
    # ω[σ,σ+1] = 1

    ω = zeros(σ, σ+1)
    g = getCoefficientsGLRK(σ)
    ω[σ,σ+1] = 1
    for i in 1:σ-1
        for j in 1:σ
            ω[i,j] = g.b[j] * g.c[j]^(i-1)
        end
    end

    TableauVSPARKsecondary(Symbol("VSPARKGLRKLobIIIAC"), o, s, σ, coeff_q, coeff_p, coeff_q̃, coeff_p̃, ω)#, get_lobatto_d_vector(σ))
end


function getTableauVSPARKGLRKLobIIIAD(s, σ=s+1)
    o = 2s

    glrk_q = getCoefficientsGLRK(s)
    glrk_p = get_symplectic_conjugate_coefficients(glrk_q)

    a_q_1 = glrk_q.a
    b_q_1 = glrk_q.b

    a_p_1 = glrk_p.a
    b_p_1 = glrk_p.b

    a_p_2 = glrk_p.a
    b_p_2 = glrk_p.b

    lobIIIA = getCoefficientsLobIIIA(σ)
    lobIIID = getCoefficientsLobIIID(σ)

    # α_q_2 and α_p_3 need to be conjugate symplectic
    # α_q_1 and α_p_2 are free

    α_q_2 = lobIIID.a
    α_p_3 = lobIIID.a
    α_q_1 = lobIIID.a
    α_p_2 = lobIIID.a

    b_q_2 = lobIIID.b
    b_p_3 = lobIIID.b

    β_q_1 = glrk_q.b
    β_q_2 = lobIIID.b

    β_p_1 = glrk_p.b
    β_p_2 = lobIIID.b
    β_p_3 = lobIIID.b

    # a_p_3_ij = b_p_3_j / b_q_1_i * (b_q_1_i - α_q_1_ji)
    a_p_3 = zeros(s,σ)
    for i in 1:s
        for j in 1:σ
            a_p_3[i,j] = b_p_3[j] / b_q_1[i] * (b_q_1[i] - α_q_1[j,i])
        end
    end

    # a_q_2_ij = b_q_2_j / b_p_2_i * (b_p_2_i - α_p_2_ji)
    a_q_2 = zeros(s,σ)
    for i in 1:s
        for j in 1:σ
            a_q_2[i,j] = b_q_2[j] / b_p_2[i] * (b_p_2[i] - α_p_2[j,i])
        end
    end

    # α_p_1_ij = b_p_1_j / b_q_2_i * (b_q_2_i - a_q_2_ji)
    α_p_1 = zeros(σ,s)
    for i in 1:σ
        for j in 1:s
            α_p_1[i,j] = b_p_1[j] / b_q_2[i] * (b_q_2[i] - a_q_2[j,i])
        end
    end

    a_q = (a_q_1, a_q_2)
    b_q = (b_q_1, b_q_2)
    c_q = glrk_q.c

    a_p = (a_p_1, a_p_2, a_p_3)
    b_p = (b_p_1, b_p_2, b_p_3)
    c_p = glrk_p.c

    α_q = (α_q_1, α_q_2)
    β_q = (β_q_1, β_q_2)
    γ_q = lobIIID.c

    α_p = (α_p_1, α_p_2, α_p_3)
    β_p = (β_p_1, β_p_2, β_p_3)
    γ_p = lobIIID.c

    coeff_q = CoefficientsSPARK(Symbol("VSPARKGLRK"),   o, s, σ, a_q, b_q, c_q)
    coeff_p = CoefficientsSPARK(Symbol("VSPARKGLRK"),   o, s, σ, a_p, b_p, c_p)
    coeff_q̃ = CoefficientsSPARK(Symbol("VSPARKLobIII"), o, σ, s, α_q, β_q, γ_q)
    coeff_p̃ = CoefficientsSPARK(Symbol("VSPARKLobIII"), o, σ, s, α_p, β_p, γ_p)

    ω = zeros(σ, σ+1)
    ω[1:σ-1,1:σ] .= lobIIIA.a[2:σ,1:σ]
    ω[σ,σ+1] = 1

    # ω = zeros(σ, σ+1)
    # g = getCoefficientsGLRK(σ)
    # ω[σ,σ+1] = 1
    # for i in 1:σ-1
    #     for j in 1:σ
    #         ω[i,j] = g.b[j] * g.c[j]^(i-1)
    #     end
    # end

    TableauVSPARKsecondary(Symbol("VSPARKGLRKLobIIIAD"), o, s, σ, coeff_q, coeff_p, coeff_q̃, coeff_p̃, ω)#, get_lobatto_d_vector(σ))
end


function getTableauVSPARKGLRKLobIIIAE(s, σ=s+1)
    o = 2s

    glrk_q = getCoefficientsGLRK(s)
    glrk_p = get_symplectic_conjugate_coefficients(glrk_q)

    a_q_1 = glrk_q.a
    b_q_1 = glrk_q.b

    a_p_1 = glrk_p.a
    b_p_1 = glrk_p.b

    a_p_2 = glrk_p.a
    b_p_2 = glrk_p.b

    lobIIIA = getCoefficientsLobIIIA(σ)
    lobIIIE = getCoefficientsLobIIIE(σ)

    # α_q_2 and α_p_3 need to be conjugate symplectic
    # α_q_1 and α_p_2 are free

    α_q_2 = lobIIIE.a
    α_p_3 = lobIIIE.a
    α_q_1 = lobIIIE.a
    α_p_2 = lobIIIE.a

    b_q_2 = lobIIIE.b
    b_p_3 = lobIIIE.b

    β_q_1 = glrk_q.b
    β_q_2 = lobIIIE.b

    β_p_1 = glrk_p.b
    β_p_2 = lobIIIE.b
    β_p_3 = lobIIIE.b

    # a_p_3_ij = b_p_3_j / b_q_1_i * (b_q_1_i - α_q_1_ji)
    a_p_3 = zeros(s,σ)
    for i in 1:s
        for j in 1:σ
            a_p_3[i,j] = b_p_3[j] / b_q_1[i] * (b_q_1[i] - α_q_1[j,i])
        end
    end

    # a_q_2_ij = b_q_2_j / b_p_2_i * (b_p_2_i - α_p_2_ji)
    a_q_2 = zeros(s,σ)
    for i in 1:s
        for j in 1:σ
            a_q_2[i,j] = b_q_2[j] / b_p_2[i] * (b_p_2[i] - α_p_2[j,i])
        end
    end

    # α_p_1_ij = b_p_1_j / b_q_2_i * (b_q_2_i - a_q_2_ji)
    α_p_1 = zeros(σ,s)
    for i in 1:σ
        for j in 1:s
            α_p_1[i,j] = b_p_1[j] / b_q_2[i] * (b_q_2[i] - a_q_2[j,i])
        end
    end

    a_q = (a_q_1, a_q_2)
    b_q = (b_q_1, b_q_2)
    c_q = glrk_q.c

    a_p = (a_p_1, a_p_2, a_p_3)
    b_p = (b_p_1, b_p_2, b_p_3)
    c_p = glrk_p.c

    α_q = (α_q_1, α_q_2)
    β_q = (β_q_1, β_q_2)
    γ_q = lobIIIE.c

    α_p = (α_p_1, α_p_2, α_p_3)
    β_p = (β_p_1, β_p_2, β_p_3)
    γ_p = lobIIIE.c

    coeff_q = CoefficientsSPARK(Symbol("VSPARKGLRK"),   o, s, σ, a_q, b_q, c_q)
    coeff_p = CoefficientsSPARK(Symbol("VSPARKGLRK"),   o, s, σ, a_p, b_p, c_p)
    coeff_q̃ = CoefficientsSPARK(Symbol("VSPARKLobIII"), o, σ, s, α_q, β_q, γ_q)
    coeff_p̃ = CoefficientsSPARK(Symbol("VSPARKLobIII"), o, σ, s, α_p, β_p, γ_p)

    ω = zeros(σ, σ+1)
    ω[1:σ-1,1:σ] .= lobIIIA.a[2:σ,1:σ]
    ω[σ,σ+1] = 1

    # ω = zeros(σ, σ+1)
    # g = getCoefficientsGLRK(σ)
    # ω[σ,σ+1] = 1
    # for i in 1:σ-1
    #     for j in 1:σ
    #         ω[i,j] = g.b[j] * g.c[j]^(i-1)
    #     end
    # end

    TableauVSPARKsecondary(Symbol("VSPARKGLRKLobIIIAE"), o, s, σ, coeff_q, coeff_p, coeff_q̃, coeff_p̃, ω)#, get_lobatto_d_vector(σ))
end
