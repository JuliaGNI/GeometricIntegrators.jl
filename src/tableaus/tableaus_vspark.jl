

function getTableauVSPARKGLRK(s)
    g = getCoefficientsGLRK(s)

    ω = zeros(s-1, s)

    for i in 1:s-1
        for j in 1:s
            ω[i,j] = g.b[j] * g.c[j]^(i-1)
        end
    end

    return TableauVSPARK(Symbol("vsparkglrk", s), g.o,
                        g.a, g.a, g.a, g.a,
                        g.a, g.a, g.a, g.a,
                        g.b, g.b, g.b, g.b,
                        g.c, g.c, g.c, g.c,
                        ω)
end


function getTableauMidpointProjection(name, q::CoefficientsRK{T}, p::CoefficientsRK{T}, d=[]; R∞=1) where {T}
    @assert q.s == p.s
    s = q.s

    g = getCoefficientsGLRK(1)
    α̃ = g.a
    β = g.b
    γ = g.c

    α  = zeros(T, s, 1)
    α .= 0.5

    q_ã = zeros(T, 1, s)
    for i in 1:s
        q_ã[1,i] = q.b[i] / β[1] * ( β[1] - α[i,1] )
    end

    p_ã = zeros(T, 1, s)
    for i in 1:s
        p_ã[1,i] = p.b[i] / β[1] * ( β[1] - α[i,1] )
    end

    β  = zero(g.b)
    β .= α̃[1,1] * (1 + R∞)

    δ = zeros(T, 1)
    ω = zeros(T, 0, 1)

    return TableauVSPARK(name, min(q.o, p.o),
                        q.a, p.a, α, α,
                        q_ã, p_ã, α̃, α̃,
                        q.b, p.b, β, β,
                        q.c, p.c, γ, δ,
                        ω)
end

"Tableau for Gauss-Legendre method with s stages and midpoint projection."
function getTableauGLRKpMidpoint(s)
    glrk = getCoefficientsGLRK(s)
    R∞ = (-1)^s
    getTableauMidpointProjection(Symbol("vsparkglrk", s, "pMidpoint"), glrk, glrk; R∞=R∞)
end



function getTableauSymmetricProjection(name, q::CoefficientsRK{T}, p::CoefficientsRK{T}, d=[]; R∞=1) where {T}

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
    c_λ = [ 0.0, 1.0]
    d_λ = [-1.0, 1.0]

    ω_λ = [0.5, R∞*0.5]
    ω_λ = reshape(ω_λ, 1, length(ω_λ))


    if length(d) == 0
        return TableauVSPARK(name, o,
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ)
    else
        @assert length(d) == q.s == p.s

        return TableauVSPARK(name, o,
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ, d)
    end

end



"Tableau for Gauss-Lobatto IIIA-IIIB method with two stages and symmetric projection."
function getTableauLobIIIAIIIB2pSymmetric()
    d = [+1.0, -1.0]

    getTableauSymmetricProjection(:LobIIIAIIIB2pSymmetric, getCoefficientsLobIIIA2(), getCoefficientsLobIIIB2(), d; R∞=-1)
end

"Tableau for Gauss-Lobatto IIIA-IIIB method with three stages and symmetric projection."
function getTableauLobIIIAIIIB3pSymmetric()
    d = [+0.5, -1.0, +0.5]

    getTableauSymmetricProjection(:LobIIIAIIIB3pSymmetric, getCoefficientsLobIIIA3(), getCoefficientsLobIIIB3(), d; R∞=+1)
end

"Tableau for Gauss-Legendre method with s stages and symplectic projection."
function getTableauGLRKpSymmetric(s)
    glrk = getCoefficientsGLRK(s)
    R∞ = (-1)^s

    getTableauSymmetricProjection(Symbol("vpglrk", s, "pSymmetric"), glrk, glrk; R∞=R∞)
end
