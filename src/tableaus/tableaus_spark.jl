
"SPARK tableau for Gauss-Legendre Runge-Kutta method with s stages."
function getTableauSPARKGLRK(s)
    g = getCoefficientsGLRK(s)
    δ = zeros(0, s)

    return TableauSPARK(Symbol("sparkglrk", s), g.o,
                        g.a, g.a, g.a, g.a,
                        g.a, g.a, g.a, g.a,
                        g.b, g.b, g.b, g.b,
                        g.c, g.c, g.c, g.b,
                        get_GLRK_ω_matrix(s), δ)
end

"SPARK tableau for Gauss-Lobatto methods."
function getTableauSPARKLob(A, B)
    @assert A.s == B.s

    s = A.s
    o = min(A.o, B.o)
    δ = zeros(0, s)

    return TableauSPARK(Symbol("sparklob", s), o,
                        A.a, B.a, A.a, B.a,
                        A.a, B.a, A.a, B.a,
                        A.b, B.b, A.b, B.b,
                        A.c, B.c, A.c, A.b,
                        get_lobatto_ω_matrix(s), δ,
                        get_lobatto_d_vector(s))
end

"SPARK tableau for Gauss-Lobatto IIIA-IIIB method with s stages."
function getTableauSPARKLobIIIAIIIB(s)
    getTableauSPARKLob(getCoefficientsLobIIIA(s), getCoefficientsLobIIIB(s))
end

"SPARK tableau for Gauss-Legendre/Gauss-Lobatto methods."
function getTableauSPARKGLRKLobIIIAIIIB(s, σ=s+1)
    g = getCoefficientsGLRK(s)
    A = getCoefficientsLobIIIA(σ)
    B = getCoefficientsLobIIIB(σ)

    α = A.a[2:σ, 1:σ]
    ã = B.a[1:σ, 1:s]

    o = min(g.o, A.o, B.o)
    δ = zeros(0, σ)

    return TableauSPARK(Symbol("spark_glrk_lob", s), o,
                        g.a, g.a, α,   α,
                        ã,   ã,   A.a, B.a,
                        g.b, g.b, A.b, B.b,
                        g.c, g.c, A.c, A.b,
                        get_lobatto_ω_matrix(σ), δ)
end

