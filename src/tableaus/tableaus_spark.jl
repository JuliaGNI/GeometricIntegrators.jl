
function getTableauSPARKGLRK(s)
    g = getCoefficientsGLRK(s)
    δ = zeros(0, s)

    return TableauSPARK(Symbol("sparkglrk", s), g.o,
                        g.a, g.a, g.a, g.a,
                        g.a, g.a, g.a, g.a,
                        g.b, g.b, g.b, g.b,
                        g.c, g.c, g.c, g.c,
                        get_GLRK_ω_matrix(s), δ)
end
