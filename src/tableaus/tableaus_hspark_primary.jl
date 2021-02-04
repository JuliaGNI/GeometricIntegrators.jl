
function TableauHSPARKSymmetricProjection(name, q::Tableau{T}, p::Tableau{T}, d=[]; R∞=1) where {T}

    @assert q.s == p.s

    o = min(q.o, p.o)

    a_q = q.a
    a_p = p.a

    α_q = zeros(T, q.s, 2)
    α_q[:,1] .= 0.5

    α_p = zeros(T, p.s, 2)
    α_p[:,1] .= 0.5

    a_q̃ = vcat(zero(q.b)', q.b')
    a_p̃ = vcat(zero(p.b)', p.b')

    α_q̃ = [[0.0  0.0]
           [0.5  0.5]]
    α_p̃ = [[0.0  0.0]
           [0.5  0.5]]

    b_q = q.b
    b_p = p.b
    β_q = [0.5, 0.5]
    β_p = [0.5, 0.5]

    c_q = q.c
    c_p = p.c
    c_λ = [ 0.0, 1.0]
    d_λ = [ 0.5, 0.5*R∞]

    ω_λ = reshape(T[1  1  0], (1,3))
    δ_λ = reshape(T[-1   R∞], (1,2))


    if length(d) == 0
        return TableauHSPARKprimary(name, o,
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ, δ_λ)
    else
        @assert length(d) == q.s == p.s

        return TableauHSPARKprimary(name, o,
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ, δ_λ, d)
    end

end

"Tableau for Gauss-Lobatto IIIA-IIIB method with s stages and symmetric projection."
function TableauHSPARKLobattoIIIAIIIBpSymmetric(s)
    TableauHSPARKSymmetricProjection(Symbol("HSPARKLobattoIIIAIIIB($s)pSymmetric"), TableauLobattoIIIA(s), TableauLobattoIIIB(s), get_lobatto_nullvector(s))
end

"Tableau for Gauss-Lobatto IIIB-IIIA method with s stages and symmetric projection."
function TableauHSPARKLobattoIIIBIIIApSymmetric(s)
    TableauHSPARKSymmetricProjection(Symbol("HSPARKLobattoIIIBIIIA($s)pSymmetric"), TableauLobattoIIIB(s), TableauLobattoIIIA(s), get_lobatto_nullvector(s))
end

"Tableau for Gauss-Legendre method with s stages and symplectic projection."
function TableauHSPARKGLRKpSymmetric(s)
    glrk = TableauGauss(s)
    TableauHSPARKSymmetricProjection(Symbol("HSPARKGLRK($s)"), glrk, glrk)
end
