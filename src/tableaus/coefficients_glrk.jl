
function get_GLRK_ω_matrix(s)
    ω = zeros(s, s+1)
    g = TableauGauss(s)
    ω[s,s+1] = 1
    for i in 1:s-1
        for j in 1:s
            ω[i,j] = g.b[j] * g.c[j]^(i-1)
        end
    end
    return ω
end
