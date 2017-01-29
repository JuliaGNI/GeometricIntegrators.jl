
function getTableauSymmetricSymplecticProjection{T}(name, q::CoefficientsRK{T}, p::CoefficientsRK{T}, d=[])

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

    ω_λ = [0.5, 0.5]
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


function getTableauSymmetricConjugateProjection{T}(name, q::CoefficientsRK{T}, p::CoefficientsRK{T}, d=[])

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
           [0.5 -0.5]]
    α_p̃ = [[0.0  0.0]
           [0.5 -0.5]]

    b_q = q.b
    b_p = p.b
    β_q = [0.5, -0.5]
    β_p = [0.5, -0.5]

    c_q = q.c
    c_p = p.c
    c_λ = [0.0, 1.0]
    d_λ = [0.0, 1.0]

    ω_λ = [0.5, 0.5]
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


function getTableauLobIIIAB2sp()
    d = [+1.0, -1.0]

    getTableauSymmetricConjugateProjection(:LobIIIAB2p, getCoefficientsLobIIIA(), getCoefficientsLobIIIB(), d)
end
