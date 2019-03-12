
function getTableauSymplecticProjection(name, q::CoefficientsRK{T}, p::CoefficientsRK{T}, d=[]; R∞=1) where {T}

    @assert q.s == p.s

    o = min(q.o, p.o)

    a_q = q.a
    a_p = p.a

    α_q = zeros(T, q.s, 2)
    α_q[:,1] .= 0.5

    α_p = zeros(T, p.s, 2)
    α_p[:,1] .= 0.5

    a_q̃ = transpose(hcat(zeros(q.b), q.b))
    a_p̃ = transpose(hcat(zeros(p.b), p.b))

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


    if length(d) == 0
        return TableauVPARK(name, o,
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ)
    else
        @assert length(d) == q.s == p.s

        return TableauVPARK(name, o,
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            d)
    end

end


"Tableau for Gauss-Lobatto IIIA-IIIB method with two stages and symplectic projection."
function getTableauLobIIIAIIIB2pSymplectic()
    d = [+1.0, -1.0]

    getTableauSymplecticProjection(:LobIIIAIIIB2pSymplectic, getCoefficientsLobIIIA2(), getCoefficientsLobIIIB2(), d; R∞=-1)
end

"Tableau for Gauss-Lobatto IIIA-IIIB method with three stages and symplectic projection."
function getTableauLobIIIAIIIB3pSymplectic()
    d = [+0.5, -1.0, +0.5]

    getTableauSymplecticProjection(:LobIIIAIIIB3pSymplectic, getCoefficientsLobIIIA3(), getCoefficientsLobIIIB3(), d; R∞=+1)
end

"Tableau for Gauss-Legendre method with s stages and symplectic projection."
function getTableauGLRKpSymplectic(s)
    glrk = getCoefficientsGLRK(s)
    R∞ = -1^s

    getTableauSymplecticProjection(Symbol("vpglrk", s, "pSymplectic"), glrk, glrk; R∞=R∞)
end
