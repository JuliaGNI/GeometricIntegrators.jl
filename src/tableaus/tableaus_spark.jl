
function getTableauSPARKGLRK(s)
    g = getCoefficientsGLRK(s)

    ω = zeros(s, s+1)
    δ = zeros(0, s)

    for i in 1:s-1
        for j in 1:s
            ω[i,j] = g.b[j] * g.c[j]^(i-1)
        end
    end

    ω[s,s+1] = 1

    return TableauSPARK(Symbol("sparkglrk", s), g.o,
                        g.a, g.a, g.a, g.a,
                        g.a, g.a, g.a, g.a,
                        g.b, g.b, g.b, g.b,
                        g.c, g.c, g.c, g.c,
                        ω, δ)
end
