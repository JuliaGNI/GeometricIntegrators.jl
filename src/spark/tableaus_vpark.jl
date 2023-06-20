
function TableauSymplecticProjection(name, q::Tableau{T}, p::Tableau{T}, d=nothing; R∞=1) where {T}

    @assert q.s == p.s

    o = min(q.o, p.o)

    a_q = q.a
    a_p = p.a

    α_q = zeros(T, q.s, 2)
    α_q[:,1] .= 0.5

    α_p = zeros(T, p.s, 2)
    α_p[:,1] .= 0.5

    a_q̃ = Array(transpose(hcat(zero(q.b), q.b)))
    a_p̃ = Array(transpose(hcat(zero(p.b), p.b)))

    α_q̃ = [[0.0  0.0]
           [0.5  R∞*0.5]]
    α_p̃ = [[0.0  0.0]
           [0.5  R∞*0.5]]

    b_q = q.b
    b_p = p.b
    β_q = [0.5, R∞*0.5]
    β_p = [0.5, R∞*0.5]

    c_q = q.c
    c_p = p.c
    c_λ = [0.0, 1.0]
    d_λ = [0.0, 1.0]


    return VPARK(TableauVPARK(name, o,
                        a_q, a_p, α_q, α_p,
                        a_q̃, a_p̃, α_q̃, α_p̃,
                        b_q, b_p, β_q, β_p,
                        c_q, c_p, c_λ, d_λ,
                        d))
end


"Tableau for Gauss-Lobatto IIIA-IIIB method with s stages and symplectic projection."
function TableauLobattoIIIAIIIBpSymplectic(s)
    TableauSymplecticProjection(Symbol("LobattoIIIAIIIB($s)pSymplectic"), TableauLobattoIIIA(s), TableauLobattoIIIB(s), get_lobatto_nullvector(s); R∞=-1^(s+1))
end

"Tableau for Gauss-Lobatto IIIB-IIIA method with s stages and symplectic projection."
function TableauLobattoIIIBIIIApSymplectic(s)
    TableauSymplecticProjection(Symbol("LobattoIIIBIIIA($s)pSymplectic"), TableauLobattoIIIB(s), TableauLobattoIIIA(s), get_lobatto_nullvector(s); R∞=-1^(s+1))
end

"Tableau for Gauss-Legendre method with s stages and symplectic projection."
function TableauGausspSymplectic(s)
    glrk = TableauGauss(s)
    TableauSymplecticProjection(Symbol("GLRK($s)pSymplectic"), glrk, glrk; R∞=(-1)^s)
end
