
function getTableauSymplecticProjection(name, q, p, d=[])

    @assert q.s == p.s

    o = min(q.o, p.o)

    a_q = q.a
    a_p = p.a

    α_q = [[0.5  0.0]
           [0.5  0.0]]
    α_p = [[0.5  0.0]
           [0.5  0.0]]

    a_q̃ = [[0.0  0.0]
           [0.5  0.5]]
    a_p̃ = [[0.0  0.0]
           [0.5  0.5]]

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
    c_λ = [0.0, 1.0]
    d_λ = [0.0, 1.0]


    if length(d) == 0
        return TableauVPARK(name, o,
                            a_q, a_p, α_q, α_p, a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ)
    else
        @assert length(d) == q.s == p.s

        return TableauVPARK(name, o,
                            a_q, a_p, α_q, α_p, a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            d)
    end

end


function getTableauLobIIIAB2p()
    d = [+1.0, -1.0]

    getTableauSymplecticProjection(:LobIIIAB2p, getCoefficientsLobIIIA(), getCoefficientsLobIIIB(), d)
end