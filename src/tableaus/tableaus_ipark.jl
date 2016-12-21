
function getTableauGLRK2symmetricProjection()

    o = 4

    glrk = getTableauGLRK2()

    a_q = glrk.a
    a_p = glrk.a

    a_q̃ = [[0.5  0.0]
           [0.5  0.0]]
    a_p̃ = [[0.5  0.0]
           [0.5  0.0]]

    α_q = [[0.0  0.0]
           [0.5  0.5]]
    α_p = [[0.0  0.0]
           [0.5  0.5]]

    α_q̃ = [[0.0  0.0]
           [0.5  0.5]]
    α_p̃ = [[0.0  0.0]
           [0.5  0.5]]

    b_q = glrk.b
    b_p = glrk.b
    β_q = [0.5, 0.5]
    β_p = [0.5, 0.5]

    c_q = glrk.c
    c_p = glrk.c
    c_λ = [0.0, 1.0]
    d_λ = [0.0, 1.0]

    return TableauIPARK(:pglrk2_psymm, o,
                        a_q, a_p, α_q, α_p, a_q̃, a_p̃, α_q̃, α_p̃,
                        b_q, b_p, β_q, β_p,
                        c_q, c_p, c_λ, d_λ)

end
