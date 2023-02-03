
"SPARK tableau for Gauss-Legendre Runge-Kutta method with s stages."
function SPARKGLRK(s; name=Symbol("GLRK($s)"))
    g = TableauGauss(s)
    ω = gauss_ω_matrix(s)
    δ = zeros(0, s)

    return SPARKMethod(TableauSPARK(name, g.o,
                        g.a, g.a, g.a, g.a,
                        g.a, g.a, g.a, g.a,
                        g.b, g.b, g.b, g.b,
                        g.c, g.c, g.c, g.b,
                        ω, δ))
end

"SPARK tableau for Gauss-Lobatto methods."
function SPARKLobatto(A, B; name=Symbol("Lobatto($s)"))
    @assert A.s == B.s

    s = A.s
    o = min(A.o, B.o)
    δ = zeros(0, s)

    return SPARKMethod(TableauSPARK(name, o,
                        A.a, B.a, A.a, B.a,
                        A.a, B.a, A.a, B.a,
                        A.b, B.b, A.b, B.b,
                        A.c, B.c, A.c, A.b,
                        lobatto_ω_matrix(s), δ,
                        get_lobatto_nullvector(s)))
end

"SPARK tableau for Gauss-Lobatto IIIA-IIIB method with s stages."
function SPARKLobattoIIIAIIIB(s)
    SPARKLobatto(TableauLobattoIIIA(s), TableauLobattoIIIĀ(s); name=Symbol("LobattoIIIAIIIB($s)"))
end

"SPARK tableau for Gauss-Lobatto IIIB-IIIA method with s stages."
function SPARKLobattoIIIBIIIA(s)
    SPARKLobatto(TableauLobattoIIIB(s), TableauLobattoIIIB̄(s); name=Symbol("LobattoIIIBIIIA($s)"))
end

"SPARK tableau for Gauss-Legendre/Gauss-Lobatto-IIIA-IIIB methods."
function SPARKGLRKLobattoIIIAIIIB(s, σ=s+1)
    g = TableauGauss(s)
    A = TableauLobattoIIIA(σ)
    B = TableauLobattoIIIB(σ)

    α = A.a[2:σ, 1:σ]
    ã = B.a[1:σ, 1:s]

    o = min(g.o, A.o, B.o)
    δ = zeros(0, σ)

    return SPARKMethod(TableauSPARK(Symbol("GLRK($s)LobattoIIIAIIIB($σ)"), o,
                        g.a, g.a, α,   α,
                        ã,   ã,   A.a, B.a,
                        g.b, g.b, A.b, B.b,
                        g.c, g.c, A.c, A.b,
                        lobatto_ω_matrix(σ), δ))
end

"SPARK tableau for Gauss-Legendre/Gauss-Lobatto-IIIB-IIIA methods."
function SPARKGLRKLobattoIIIBIIIA(s, σ=s+1)
    g = TableauGauss(s)
    A = TableauLobattoIIIA(σ)
    B = TableauLobattoIIIB(σ)

    α = B.a[2:σ, 1:σ]
    ã = A.a[1:σ, 1:s]

    o = min(g.o, B.o, A.o)
    δ = zeros(0, σ)

    return SPARKMethod(TableauSPARK(Symbol("GLRK($s)LobattoIIIAIIIB($σ)"), o,
                        g.a, g.a, α,   α,
                        ã,   ã,   B.a, A.a,
                        g.b, g.b, B.b, A.b,
                        g.c, g.c, A.c, A.b,
                        lobatto_ω_matrix(σ), δ))
end



function SPARKLobatto(name, l1::Tableau{T}, l2::Tableau{T},
                                   l3::Tableau{T}, l4::Tableau{T}=l3, d=nothing; R∞=1) where {T}

    @assert l1.s == l2.s == l3.s == l4.s

    ω = lobatto_ω_matrix(l3.s)
    δ = zeros(T, 0, 1)

    return SPARKMethod(TableauSPARK(name, min(l1.o, l2.o, l3.o, l4.o),
                        l3.a, l1.a, l4.a, l2.a,
                        l3.a, l1.a, l4.a, l2.a,
                        l3.b, l1.b, l4.b, l2.b,
                        l3.c, l1.c, l4.c, l4.b,
                        ω, δ, d))
end

"Tableau for Gauss-Lobatto IIIA-IIIB-IIIC method with s stages."
function SPARKLobABC(s)
    loba = TableauLobattoIIIA(s)
    lobb = TableauLobattoIIIB(s)
    lobc = TableauLobattoIIIC(s)
    SPARKLobatto(Symbol("LobattoIIIABC($s)"), loba, lobc, lobb; R∞=(-1)^(s+1))
end

"Tableau for Gauss-Lobatto IIIA-IIIB-IIID method with s stages."
function SPARKLobABD(s)
    loba = TableauLobattoIIIA(s)
    lobb = TableauLobattoIIIB(s)
    lobd = TableauLobattoIIID(s)
    SPARKLobatto(Symbol("LobattoIIIABD($s)"), loba, lobd, lobb; R∞=(-1)^(s+1))
end

"SPARK Tableau for Variational Partitioned Runge-Kutta Methods."
function SPARKVPRK(name, q::Tableau{T}, p::Tableau{T}, d=nothing; R∞=1) where {T}

    @assert q.s == p.s

    ω = zeros(T, q.s, q.s+1)
    δ = zeros(T, 0, 1)

    for i in 1:q.s
        ω[i,i] = 1
    end

    return SPARKMethod(TableauSPARK(name, min(q.o, p.o),
                        q.a, p.a, q.a, p.a,
                        q.a, p.a, q.a, p.a,
                        q.b, p.b, q.b, p.b,
                        q.c, p.c, q.c, q.b,
                        ω, δ, d))
end

"Tableau for Variational Gauss-Legendre method with s stages."
function SPARKGLVPRK(s)
    glrk = TableauGauss(s)
    SPARKVPRK(Symbol("GLVPRK($s)"), glrk, glrk; R∞=(-1)^s)
end
