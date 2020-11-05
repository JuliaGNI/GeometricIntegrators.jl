
function getTableauVSPARKInternalProjection(name, q::CoefficientsRK{T}, p::CoefficientsRK{T}, d=[]; R∞=1) where {T}

    @assert q.s == p.s

    s = q.s
    o = min(q.o, p.o)

    α_q = zeros(T, s, s)
    α_p = zeros(T, s, s)

    for i in 1:s
        α_q[i,:] .= q.b ./ 2
        α_p[i,:] .= p.b ./ 2
    end

    β_q = q.b .* (1 + R∞) ./ 2
    β_p = p.b .* (1 + R∞) ./ 2
    
    c_λ = q.c
    d_λ = q.b
    ω_λ = zeros(T, 1, s+1)
    δ_λ = zeros(T, s-1, s)
    
    ω_λ[1,s+1] = 1
    
    for i in 1:s-1
        δ_λ[i,i] = +1
        δ_λ[i,s] = -1
    end


    if length(d) == 0
        return TableauVSPARKprimary(name, o,
                            q.a, p.a, α_q, α_p,
                            q.a, p.a, α_q, α_p,
                            q.b, p.b, β_q, β_p,
                            q.c, p.c, c_λ, d_λ,
                            ω_λ, δ_λ)
    else
        @assert length(d) == q.s == p.s

        return TableauVSPARKprimary(name, o,
                            q.a, p.a, α_q, α_p,
                            q.a, p.a, α_q, α_p,
                            q.b, p.b, β_q, β_p,
                            q.c, p.c, c_λ, d_λ,
                            ω_λ, δ_λ, d)
    end

end

"Tableau for Gauss-Legendre method with s stages and symplectic projection."
function getTableauVSPARKGLRKpInternal(s)
    glrk = getCoefficientsGLRK(s)
    getTableauVSPARKInternalProjection(Symbol("vpglrk", s, "pInternal"), glrk, glrk; R∞=(-1)^s)
end




function getTableauVSPARKMidpointProjection(name, q::CoefficientsRK{T}, p::CoefficientsRK{T}, d=[]; R∞=1) where {T}
    @assert q.s == p.s
    s = q.s

    g = getCoefficientsGLRK(1)
    α̃ = g.a
    β = g.b
    γ = g.c

    α = 0.5 * ones(T, s, 1)

    q_ã = zeros(T, 1, s)
    for i in 1:s
        q_ã[1,i] = q.b[i] / β[1] * ( β[1] - α[i,1] )
    end

    p_ã = zeros(T, 1, s)
    for i in 1:s
        p_ã[1,i] = p.b[i] / β[1] * ( β[1] - α[i,1] )
    end

    β .= (1 + R∞) / 2
    d  = ones(T, 1)
    ω  = reshape(T[0  1], (1,2))
    δ  = zeros(T, 0, 1)

    return TableauVSPARKprimary(name, min(q.o, p.o),
                        q.a, p.a, α, α,
                        q_ã, p_ã, α̃, α̃,
                        q.b, p.b, β, β,
                        q.c, p.c, γ, d,
                        ω, δ)
end

"Tableau for Gauss-Legendre method with s stages and midpoint projection."
function getTableauVSPARKGLRKpMidpoint(s)
    glrk = getCoefficientsGLRK(s)
    getTableauVSPARKMidpointProjection(Symbol("vsparkglrk", s, "pMidpoint"), glrk, glrk; R∞=(-1)^s)
end



function getTableauVSPARKSymmetricProjection(name, q::CoefficientsRK{T}, p::CoefficientsRK{T}, d=[]; R∞=1) where {T}

    @assert q.s == p.s

    o = min(q.o, p.o)

    a_q = q.a
    a_p = p.a

    α_q = zeros(T, q.s, 2)
    α_q[:,1] .= 0.5

    α_p = zeros(T, p.s, 2)
    α_p[:,1] .= 0.5

    a_q̃ = vcat(zero(q.b)', q.b')
    a_p̃ = vcat(zero(p.b)', p.b')

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
    d_λ = [ 0.5, 0.5]

    ω_λ = reshape(T[1  R∞  0], (1,3))
    δ_λ = reshape(T[-1  +1], (1,2))


    if length(d) == 0
        return TableauVSPARKprimary(name, o,
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ, δ_λ)
    else
        @assert length(d) == q.s == p.s

        return TableauVSPARKprimary(name, o,
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ, δ_λ, d)
    end

end


"Tableau for Gauss-Lobatto IIIA-IIIB method with two stages and symmetric projection."
function getTableauVSPARKLobIIIAIIIB2pSymmetric()
    getTableauVSPARKSymmetricProjection(:LobIIIAIIIB2pSymmetric, getCoefficientsLobIIIA2(), getCoefficientsLobIIIB2(), get_lobatto_d_vector(2); R∞=-1)
end

"Tableau for Gauss-Lobatto IIIA-IIIB method with three stages and symmetric projection."
function getTableauVSPARKLobIIIAIIIB3pSymmetric()
    getTableauVSPARKSymmetricProjection(:LobIIIAIIIB3pSymmetric, getCoefficientsLobIIIA3(), getCoefficientsLobIIIB3(), get_lobatto_d_vector(3); R∞=+1)
end

"Tableau for Gauss-Lobatto IIIA-IIIB method with four stages and symmetric projection."
function getTableauVSPARKLobIIIAIIIB4pSymmetric()
    getTableauVSPARKSymmetricProjection(:LobIIIAIIIB4pSymmetric, getCoefficientsLobIIIA4(), getCoefficientsLobIIIB4(), get_lobatto_d_vector(4); R∞=-1)
end

"Tableau for Gauss-Lobatto IIIA-IIIB method with five stages and symmetric projection."
function getTableauVSPARKLobIIIAIIIB5pSymmetric()
    getTableauVSPARKSymmetricProjection(:LobIIIAIIIB5pSymmetric, getCoefficientsLobIIIA5(), getCoefficientsLobIIIB5(), get_lobatto_d_vector(5); R∞=+1)
end

"Tableau for Gauss-Legendre method with s stages and symplectic projection."
function getTableauVSPARKGLRKpSymmetric(s)
    glrk = getCoefficientsGLRK(s)
    getTableauVSPARKSymmetricProjection(Symbol("vpglrk", s, "pSymmetric"), glrk, glrk; R∞=(-1)^s)
end



function getTableauVSPARKSymmetricLobABProjection(name, q::CoefficientsRK{T}, p::CoefficientsRK{T}, d=[]; R∞=1) where {T}

    @assert q.s == p.s

    o = min(q.o, p.o)

    loba = getCoefficientsLobIIIA2()
    lobb = getCoefficientsLobIIIB2()

    a_q = q.a
    a_p = p.a
    b_q = q.b
    b_p = p.b

    a_q̃ = vcat(zero(q.b)', q.b')
    a_p̃ = vcat(zero(p.b)', p.b')
#    a_p̃ = vcat(p.b' ./ 2, p.b' ./ 2)

    β_q = [0.5, R∞*0.5]
    β_p = [0.5, R∞*0.5]

    α_q̃ = loba.a
    α_p̃ = lobb.a

    α_q = zeros(T, q.s, 2)
    α_p = zeros(T, p.s, 2)

    for i in 1:q.s
        for j in 1:2
            α_q[i,j] = β_q[j] / b_p[i] * ( b_p[i] - a_p̃[j,i] )
            α_p[i,j] = β_p[j] / b_q[i] * ( b_q[i] - a_q̃[j,i] )
        end
    end

    c_q = q.c
    c_p = p.c
    c_λ = [ 0.0, 1.0]
    d_λ = [ 0.5, 0.5]

    ω_λ = reshape(T[0  0  1], (1,3))
    δ_λ = reshape(T[-1  +1], (1,2))


    if length(d) == 0
        return TableauVSPARKprimary(name, o,
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ, δ_λ)
    else
        @assert length(d) == q.s == p.s

        return TableauVSPARKprimary(name, o,
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ, δ_λ, d)
    end

end


"Tableau for Gauss-Lobatto IIIA-IIIB method with two stages and symmetric projection."
function getTableauVSPARKLobIIIAIIIB2pSymmetricLobAB()
    getTableauVSPARKSymmetricLobABProjection(:LobIIIAIIIB2pSymmetricLobAB, getCoefficientsLobIIIA2(), getCoefficientsLobIIIB2(), get_lobatto_d_vector(2); R∞=-1)
end

"Tableau for Gauss-Lobatto IIIA-IIIB method with three stages and symmetric projection."
function getTableauVSPARKLobIIIAIIIB3pSymmetricLobAB()
    getTableauVSPARKSymmetricLobABProjection(:LobIIIAIIIB3pSymmetricLobAB, getCoefficientsLobIIIA3(), getCoefficientsLobIIIB3(), get_lobatto_d_vector(3); R∞=+1)
end

"Tableau for Gauss-Lobatto IIIA-IIIB method with four stages and symmetric projection."
function getTableauVSPARKLobIIIAIIIB4pSymmetricLobAB()
    getTableauVSPARKSymmetricLobABProjection(:LobIIIAIIIB4pSymmetricLobAB, getCoefficientsLobIIIA4(), getCoefficientsLobIIIB4(), get_lobatto_d_vector(4); R∞=-1)
end

"Tableau for Gauss-Lobatto IIIA-IIIB method with five stages and symmetric projection."
function getTableauVSPARKLobIIIAIIIB5pSymmetricLobAB()
    getTableauVSPARKSymmetricLobABProjection(:LobIIIAIIIB5pSymmetricLobAB, getCoefficientsLobIIIA5(), getCoefficientsLobIIIB5(), get_lobatto_d_vector(5); R∞=+1)
end

"Tableau for Gauss-Legendre method with s stages and symplectic projection."
function getTableauVSPARKGLRKpSymmetricLobAB(s)
    glrk = getCoefficientsGLRK(s)
    getTableauVSPARKSymmetricLobABProjection(Symbol("vpglrk", s, "pSymmetricLobAB"), glrk, glrk; R∞=(-1)^s)
end



function getTableauVSPARKSymmetricLobBAProjection(name, q::CoefficientsRK{T}, p::CoefficientsRK{T}, d=[]; R∞=1) where {T}

    @assert q.s == p.s

    o = min(q.o, p.o)

    loba = getCoefficientsLobIIIA2()
    lobb = getCoefficientsLobIIIB2()

    a_q = q.a
    a_p = p.a
    b_q = q.b
    b_p = p.b

    a_q̃ = vcat(zero(q.b)', q.b')
    a_p̃ = vcat(zero(p.b)', p.b')
#    a_p̃ = vcat(p.b' ./ 2, p.b' ./ 2)

    β_q = [0.5, R∞*0.5]
    β_p = [0.5, R∞*0.5]

    α_q̃ = lobb.a
    α_p̃ = loba.a

    α_q = zeros(T, q.s, 2)
    α_p = zeros(T, p.s, 2)

    for i in 1:q.s
        for j in 1:2
            α_q[i,j] = β_q[j] / b_p[i] * ( b_p[i] - a_p̃[j,i] )
            α_p[i,j] = β_p[j] / b_q[i] * ( b_q[i] - a_q̃[j,i] )
        end
    end

    c_q = q.c
    c_p = p.c
    c_λ = [ 0.0, 1.0]
    d_λ = [ 0.5, 0.5]

    ω_λ = reshape(T[0  0  1], (1,3))
    δ_λ = reshape(T[-1  +1], (1,2))


    if length(d) == 0
        return TableauVSPARKprimary(name, o,
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ, δ_λ)
    else
        @assert length(d) == q.s == p.s

        return TableauVSPARKprimary(name, o,
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ, δ_λ, d)
    end

end


"Tableau for Gauss-Lobatto IIIA-IIIB method with two stages and symmetric projection."
function getTableauVSPARKLobIIIAIIIB2pSymmetricLobBA()
    getTableauVSPARKSymmetricLobBAProjection(:LobIIIAIIIB2pSymmetricLobBA, getCoefficientsLobIIIA2(), getCoefficientsLobIIIB2(), get_lobatto_d_vector(2); R∞=-1)
end

"Tableau for Gauss-Lobatto IIIA-IIIB method with three stages and symmetric projection."
function getTableauVSPARKLobIIIAIIIB3pSymmetricLobBA()
    getTableauVSPARKSymmetricLobBAProjection(:LobIIIAIIIB3pSymmetricLobBA, getCoefficientsLobIIIA3(), getCoefficientsLobIIIB3(), get_lobatto_d_vector(3); R∞=+1)
end

"Tableau for Gauss-Lobatto IIIA-IIIB method with four stages and symmetric projection."
function getTableauVSPARKLobIIIAIIIB4pSymmetricLobBA()
    getTableauVSPARKSymmetricLobBAProjection(:LobIIIAIIIB4pSymmetricLobBA, getCoefficientsLobIIIA4(), getCoefficientsLobIIIB4(), get_lobatto_d_vector(4); R∞=-1)
end

"Tableau for Gauss-Lobatto IIIA-IIIB method with five stages and symmetric projection."
function getTableauVSPARKLobIIIAIIIB5pSymmetricLobBA()
    getTableauVSPARKSymmetricLobBAProjection(:LobIIIAIIIB5pSymmetricLobBA, getCoefficientsLobIIIA5(), getCoefficientsLobIIIB5(), get_lobatto_d_vector(5); R∞=+1)
end

"Tableau for Gauss-Legendre method with s stages and symplectic projection."
function getTableauVSPARKGLRKpSymmetricLobBA(s)
    glrk = getCoefficientsGLRK(s)
    getTableauVSPARKSymmetricLobBAProjection(Symbol("vpglrk", s, "pSymmetricLobBA"), glrk, glrk; R∞=(-1)^s)
end



function getTableauVSPARKLobIIIAIIIBProjection(name, q::CoefficientsRK{T}, p::CoefficientsRK{T}, d=[]; R∞=1) where {T}

    @assert q.s == p.s

    s = q.s
    o = min(q.o, p.o)

    loba = getCoefficientsLobIIIA2()
    lobb = getCoefficientsLobIIIB2()

    a_q = q.a
    a_p = p.a
    b_q = q.b
    b_p = p.b

    # α_q = hcat(ones(s)/2, zeros(s))
    # α_p = hcat(ones(s)/2, zeros(s))

    α_q = zeros(T, s, 2)
    α_p = zeros(T, s, 2)
    for i in 1:s
        α_q[i,:] .= loba.b
        α_p[i,:] .= lobb.b
    end

    β_q = loba.b
    β_p = lobb.b

    α_q̃ = loba.a
    α_p̃ = lobb.a

    a_q̃ = zeros(T, 2, s)
    a_p̃ = zeros(T, 2, s)

    for i in 1:2
        for j in 1:s
            a_q̃[i,j] = b_q[j] / β_p[i] * ( β_p[i] - α_p[j,i] )
            a_p̃[i,j] = b_p[j] / β_q[i] * ( β_q[i] - α_q[j,i] )
        end
    end

    c_q = q.c
    c_p = p.c
    c_λ = loba.c
    d_λ = loba.b

    ω_λ = T[0.5  0.5  0.0
            0.0  0.0  1.0]
    δ_λ = zeros(T, 0, 1)


    if length(d) == 0
        return TableauVSPARKprimary(name, o,
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ, δ_λ)
    else
        @assert length(d) == q.s == p.s

        return TableauVSPARKprimary(name, o,
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ, δ_λ, d)
    end

end

"Tableau for Gauss-Lobatto IIIA-IIIB method with two stages and symmetric projection."
function getTableauVSPARKLobIIIAIIIB2pLobIIIAIIIB()
    getTableauVSPARKLobIIIAIIIBProjection(:LobIIIAIIIB2pLobIIIAIIIB, getCoefficientsLobIIIA2(), getCoefficientsLobIIIB2(), get_lobatto_d_vector(2); R∞=-1)
end

"Tableau for Gauss-Lobatto IIIA-IIIB method with three stages and symmetric projection."
function getTableauVSPARKLobIIIAIIIB3pLobIIIAIIIB()
    getTableauVSPARKLobIIIAIIIBProjection(:LobIIIAIIIB3pLobIIIAIIIB, getCoefficientsLobIIIA3(), getCoefficientsLobIIIB3(), get_lobatto_d_vector(3); R∞=+1)
end

"Tableau for Gauss-Lobatto IIIA-IIIB method with four stages and symmetric projection."
function getTableauVSPARKLobIIIAIIIB4pLobIIIAIIIB()
    getTableauVSPARKLobIIIAIIIBProjection(:LobIIIAIIIB4pLobIIIAIIIB, getCoefficientsLobIIIA4(), getCoefficientsLobIIIB4(), get_lobatto_d_vector(4); R∞=-1)
end

"Tableau for Gauss-Lobatto IIIA-IIIB method with five stages and symmetric projection."
function getTableauVSPARKLobIIIAIIIB5pLobIIIAIIIB()
    getTableauVSPARKLobIIIAIIIBProjection(:LobIIIAIIIB5pLobIIIAIIIB, getCoefficientsLobIIIA5(), getCoefficientsLobIIIB5(), get_lobatto_d_vector(5); R∞=+1)
end

"Tableau for Gauss-Legendre method with s stages and symplectic projection."
function getTableauVSPARKGLRKpLobIIIAIIIB(s)
    glrk = getCoefficientsGLRK(s)
    getTableauVSPARKLobIIIAIIIBProjection(Symbol("vpglrk", s, "pLobIIIAIIIB"), glrk, glrk; R∞=(-1)^s)
end




function getTableauVSPARKLobIIIBIIIAProjection(name, q::CoefficientsRK{T}, p::CoefficientsRK{T}, d=[]; R∞=1) where {T}

    @assert q.s == p.s

    s = q.s
    o = min(q.o, p.o)

    loba = getCoefficientsLobIIIA2()
    lobb = getCoefficientsLobIIIB2()

    a_q = q.a
    a_p = p.a
    b_q = q.b
    b_p = p.b

    # α_q = hcat(ones(s)/2, zeros(s))
    # α_p = hcat(ones(s)/2, zeros(s))

    α_q = zeros(T, s, 2)
    α_p = zeros(T, s, 2)
    for i in 1:s
        α_q[i,:] .= lobb.b
        α_p[i,:] .= loba.b
    end

    β_q = lobb.b
    β_p = loba.b

    α_q̃ = lobb.a
    α_p̃ = loba.a

    a_q̃ = zeros(T, 2, s)
    a_p̃ = zeros(T, 2, s)

    for i in 1:2
        for j in 1:s
            a_q̃[i,j] = b_q[j] / β_p[i] * ( β_p[i] - α_p[j,i] )
            a_p̃[i,j] = b_p[j] / β_q[i] * ( β_q[i] - α_q[j,i] )
        end
    end

    c_q = q.c
    c_p = p.c
    c_λ = loba.c
    d_λ = loba.b

    ω_λ = T[0.5  0.5  0.0
            0.0  0.0  1.0]
    δ_λ = zeros(T, 0, 1)


    if length(d) == 0
        return TableauVSPARKprimary(name, o,
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ, δ_λ)
    else
        @assert length(d) == q.s == p.s

        return TableauVSPARKprimary(name, o,
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ, δ_λ, d)
    end

end

"Tableau for Gauss-Lobatto IIIA-IIIB method with two stages and symmetric projection."
function getTableauVSPARKLobIIIAIIIB2pLobIIIBIIIA()
    getTableauVSPARKLobIIIBIIIAProjection(:LobIIIAIIIB2pLobIIIBIIIA, getCoefficientsLobIIIA2(), getCoefficientsLobIIIB2(), get_lobatto_d_vector(2); R∞=-1)
end

"Tableau for Gauss-Lobatto IIIA-IIIB method with three stages and symmetric projection."
function getTableauVSPARKLobIIIAIIIB3pLobIIIBIIIA()
    getTableauVSPARKLobIIIBIIIAProjection(:LobIIIAIIIB3pLobIIIBIIIA, getCoefficientsLobIIIA3(), getCoefficientsLobIIIB3(), get_lobatto_d_vector(3); R∞=+1)
end

"Tableau for Gauss-Lobatto IIIA-IIIB method with four stages and symmetric projection."
function getTableauVSPARKLobIIIAIIIB4pLobIIIBIIIA()
    getTableauVSPARKLobIIIBIIIAProjection(:LobIIIAIIIB4pLobIIIBIIIA, getCoefficientsLobIIIA4(), getCoefficientsLobIIIB4(), get_lobatto_d_vector(4); R∞=-1)
end

"Tableau for Gauss-Lobatto IIIA-IIIB method with five stages and symmetric projection."
function getTableauVSPARKLobIIIAIIIB5pLobIIIBIIIA()
    getTableauVSPARKLobIIIBIIIAProjection(:LobIIIAIIIB5pLobIIIBIIIA, getCoefficientsLobIIIA5(), getCoefficientsLobIIIB5(), get_lobatto_d_vector(5); R∞=+1)
end

"Tableau for Gauss-Legendre method with s stages and symplectic projection."
function getTableauVSPARKGLRKpLobIIIBIIIA(s)
    glrk = getCoefficientsGLRK(s)
    getTableauVSPARKLobIIIBIIIAProjection(Symbol("vpglrk", s, "pLobIIIBIIIA"), glrk, glrk; R∞=(-1)^s)
end




function getTableauVSPARKSymplecticProjection(name, q::CoefficientsRK{T}, p::CoefficientsRK{T}, d=[]; R∞=+1) where {T}

    la = getCoefficientsLobIIIA2()
    lb = getCoefficientsLobIIIB2()

    @assert q.s == p.s
    @assert q.b == p.b
    @assert q.c == p.c

    @assert la.s == lb.s
    @assert la.b == lb.b
#    @assert la.c == lb.c

    s = q.s
    σ = la.s
    ρ = 0

    a_q = q.a
    a_p = p.a
    b_q = q.b
    b_p = p.b

    # β_q = la.b
    # β_p = lb.b
    β_q = [0.5, R∞*0.5]
    β_p = [0.5, R∞*0.5]

    α_q = zeros(T, s, 2)
    α_q[:,1] .= 0.5

    α_p = zeros(T, s, 2)
    α_p[:,1] .= 0.5

#    α_q = zeros(T, s, σ)
#    α_p = zeros(T, s, σ)
#    for i in 1:s
#        for j in 1:σ
#            α_q[i,j] = #q.b[i] / β[1] * ( β[1] - α[i,1] )
#            α_p[i,j] = #p.b[i] / β[1] * ( β[1] - α[i,1] )
#        end
#    end

    a_q̃ = zeros(T, σ, s)
    a_p̃ = zeros(T, σ, s)
    for i in 1:σ
        for j in 1:s
            a_q̃[i,j] = b_q[j] / β_p[i] * (β_p[i] - α_p[j,i])
            a_p̃[i,j] = b_p[j] / β_q[i] * (β_q[i] - α_q[j,i])
        end
    end

    α_q̃ = la.a
    α_p̃ = lb.a

    c_q = q.c
    c_p = p.c
    c_λ = la.c
    d_λ = [0.5, 0.5]


    ω_λ = [0.5 0.5 0.0
           0.0 0.0 1.0]

    δ_λ = zeros(T, ρ, σ)


    if length(d) == 0
        return TableauVSPARKprimary(name, min(q.o, p.o),
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ, δ_λ)
    else
        @assert length(d) == s

        return TableauVSPARKprimary(name, min(q.o, p.o),
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ, δ_λ, d)
    end

end

function getTableauVSPARKGLRKpSymplectic(s)
    glrk = getCoefficientsGLRK(s)
    getTableauVSPARKSymplecticProjection(Symbol("VSPARK", s, "pSymplectic"), glrk, glrk; R∞=(-1)^s)
end


function getTableauVSPARKLobABCCD(s=2, T=Float64)
    lob  = getCoefficientsLobIII(s)
    loba = getCoefficientsLobIIIA(s)
    lobb = getCoefficientsLobIIIB(s)
    lobc = getCoefficientsLobIIIC(s)
    lobd = getCoefficientsLobIIID(s)
    
    a_q = lobd.a
    a_p = lobd.a
    
    α_q = lobb.a
    α_p = lobb.a

    a_q̃ = loba.a
    a_p̃ = loba.a
    
    α_q̃ = lob.a
    α_p̃ = lobc.a

    
    b_q = lobd.b
    b_p = lobd.b

    β_q̃ = lobd.b
    β_p̃ = lobd.b
    
    c_q = lobd.c
    c_p = lobd.c

    γ = lobd.c

    d = ones(T, s)
    ω = get_lobatto_ω_matrix(s)
    δ = zeros(T, 0, s)

    return TableauVSPARKprimary(:vspark_lobabccd, 2*s-2,
                        a_q, a_p, α_q, α_p,
                        a_q̃, a_p̃, α_q̃, α_p̃,
                        b_q, b_p, β_q̃, β_p̃,
                        c_q, c_p, γ, d,
                        ω, δ)
end

function getTableauVSPARKLobABCCE(s=2, T=Float64)
    lob  = getCoefficientsLobIII(s)
    loba = getCoefficientsLobIIIA(s)
    lobb = getCoefficientsLobIIIB(s)
    lobc = getCoefficientsLobIIIC(s)
    lobd = getCoefficientsLobIIID(s)
    lobe = getCoefficientsLobIIIE(s)
    
    a_q = lobe.a
    a_p = lobe.a
    
    α_q = lobb.a
    α_p = lobb.a

    a_q̃ = loba.a
    a_p̃ = loba.a
    
    α_q̃ = lob.a
    α_p̃ = lobc.a

    
    b_q = lobd.b
    b_p = lobd.b

    β_q̃ = lobd.b
    β_p̃ = lobd.b
    
    c_q = lobd.c
    c_p = lobd.c

    γ = lobd.c

    d = ones(T, s)
    ω = get_lobatto_ω_matrix(s)
    δ = zeros(T, 0, s)

    return TableauVSPARKprimary(:vspark_lobabccd, 2*s-2,
                        a_q, a_p, α_q, α_p,
                        a_q̃, a_p̃, α_q̃, α_p̃,
                        b_q, b_p, β_q̃, β_p̃,
                        c_q, c_p, γ, d,
                        ω, δ)
end

function getTableauVSPARKLobABED(s=2, T=Float64)
    lob  = getCoefficientsLobIII(s)
    loba = getCoefficientsLobIIIA(s)
    lobb = getCoefficientsLobIIIB(s)
    lobc = getCoefficientsLobIIIC(s)
    lobd = getCoefficientsLobIIID(s)
    lobe = getCoefficientsLobIIIE(s)

    a_q = lobe.a
    a_p = lobe.a
    
    α_q = lobb.a
    α_p = lobb.a

    a_q̃ = loba.a
    a_p̃ = loba.a
    
    α_q̃ = lobd.a
    α_p̃ = lobd.a

    
    b_q = lob.b
    b_p = lob.b

    β_q̃ = lob.b
    β_p̃ = lob.b
    
    c_q = lob.c
    c_p = lob.c

    γ   = lob.c

    d = ones(T, s)
    ω = get_lobatto_ω_matrix(s)
    δ = zeros(T, 0, s)

    return TableauVSPARKprimary(:vspark_lobabccd, 2*s-2,
                        a_q, a_p, α_q, α_p,
                        a_q̃, a_p̃, α_q̃, α_p̃,
                        b_q, b_p, β_q̃, β_p̃,
                        c_q, c_p, γ, d,
                        ω, δ)
end

function getTableauVSPARKLobABDE(s=2, T=Float64)
    lob  = getCoefficientsLobIII(s)
    loba = getCoefficientsLobIIIA(s)
    lobb = getCoefficientsLobIIIB(s)
    lobc = getCoefficientsLobIIIC(s)
    lobd = getCoefficientsLobIIID(s)
    lobe = getCoefficientsLobIIIE(s)

    a_q = lobd.a
    a_p = lobd.a
    
    α_q = lobb.a
    α_p = lobb.a

    a_q̃ = loba.a
    a_p̃ = loba.a
    
    α_q̃ = lobe.a
    α_p̃ = lobe.a

    
    b_q = lob.b
    b_p = lob.b

    β_q̃ = lob.b
    β_p̃ = lob.b
    
    c_q = lob.c
    c_p = lob.c

    γ   = lob.c

    d = ones(T, s)
    ω = get_lobatto_ω_matrix(s)
    δ = zeros(T, 0, s)

    return TableauVSPARKprimary(:vspark_lobabccd, 2*s-2,
                        a_q, a_p, α_q, α_p,
                        a_q̃, a_p̃, α_q̃, α_p̃,
                        b_q, b_p, β_q̃, β_p̃,
                        c_q, c_p, γ, d,
                        ω, δ)
end

function getTableauVSPARKLobABD(s=2, T=Float64)
    lob  = getCoefficientsLobIII(s)
    loba = getCoefficientsLobIIIA(s)
    lobb = getCoefficientsLobIIIB(s)
    lobc = getCoefficientsLobIIIC(s)
    lobd = getCoefficientsLobIIID(s)
    lobe = getCoefficientsLobIIIE(s)


    a_q = lobd.a
    a_p = lobd.a
    
    α_q = lobb.a
    α_p = lobb.a

    a_q̃ = loba.a
    a_p̃ = loba.a
    
    α_q̃ = loba.a
    α_p̃ = lobb.a

    
    b_q = lob.b
    b_p = lob.b

    β_q̃ = lob.b
    β_p̃ = lob.b
    
    c_q = lob.c
    c_p = lob.c

    γ   = lob.c

    d = ones(T, s)
    ω = get_lobatto_ω_matrix(s)
    δ = zeros(T, 0, s)

    return TableauVSPARKprimary(:vspark_lobabccd, 2*s-2,
                        a_q, a_p, α_q, α_p,
                        a_q̃, a_p̃, α_q̃, α_p̃,
                        b_q, b_p, β_q̃, β_p̃,
                        c_q, c_p, γ, d,
                        ω, δ)
end

function getTableauVSPARKLobABE(s=2, T=Float64)
    lob  = getCoefficientsLobIII(s)
    loba = getCoefficientsLobIIIA(s)
    lobb = getCoefficientsLobIIIB(s)
    lobc = getCoefficientsLobIIIC(s)
    lobd = getCoefficientsLobIIID(s)
    lobe = getCoefficientsLobIIIE(s)


    a_q = lobe.a
    a_p = lobe.a
    
    α_q = lobb.a
    α_p = lobb.a

    a_q̃ = loba.a
    a_p̃ = loba.a
    
    α_q̃ = loba.a
    α_p̃ = lobb.a

    
    b_q = lob.b
    b_p = lob.b

    β_q̃ = lob.b
    β_p̃ = lob.b
    
    c_q = lob.c
    c_p = lob.c

    γ   = lob.c

    d = ones(T, s)
    ω = get_lobatto_ω_matrix(s)
    δ = zeros(T, 0, s)

    return TableauVSPARKprimary(:vspark_lobabccd, 2*s-2,
                        a_q, a_p, α_q, α_p,
                        a_q̃, a_p̃, α_q̃, α_p̃,
                        b_q, b_p, β_q̃, β_p̃,
                        c_q, c_p, γ, d,
                        ω, δ)
end

function getTableauVSPARKLobDE(s=2, T=Float64)
    lob  = getCoefficientsLobIII(s)
    loba = getCoefficientsLobIIIA(s)
    lobb = getCoefficientsLobIIIB(s)
    lobc = getCoefficientsLobIIIC(s)
    lobd = getCoefficientsLobIIID(s)
    lobe = getCoefficientsLobIIIE(s)


    a_q = lobd.a
    a_p = lobd.a
    
    α_q = lobe.a
    α_p = lobe.a

    a_q̃ = lobe.a
    a_p̃ = lobe.a
    
    α_q̃ = lobd.a
    α_p̃ = lobd.a

    
    b_q = lob.b
    b_p = lob.b

    β_q̃ = lob.b
    β_p̃ = lob.b
    
    c_q = lob.c
    c_p = lob.c

    γ   = lob.c

    d = ones(T, s)
    ω = get_lobatto_ω_matrix(s)
    δ = zeros(T, 0, s)

    return TableauVSPARKprimary(:vspark_lobabccd, 2*s-2,
                        a_q, a_p, α_q, α_p,
                        a_q̃, a_p̃, α_q̃, α_p̃,
                        b_q, b_p, β_q̃, β_p̃,
                        c_q, c_p, γ, d,
                        ω, δ)
end

function getTableauVSPARKLobED(s=2, T=Float64)
    lob  = getCoefficientsLobIII(s)
    loba = getCoefficientsLobIIIA(s)
    lobb = getCoefficientsLobIIIB(s)
    lobc = getCoefficientsLobIIIC(s)
    lobd = getCoefficientsLobIIID(s)
    lobe = getCoefficientsLobIIIE(s)


    a_q = lobe.a
    a_p = lobe.a
    
    α_q = lobd.a
    α_p = lobd.a

    a_q̃ = lobd.a
    a_p̃ = lobd.a
    
    α_q̃ = lobe.a
    α_p̃ = lobe.a

    
    b_q = lob.b
    b_p = lob.b

    β_q̃ = lob.b
    β_p̃ = lob.b
    
    c_q = lob.c
    c_p = lob.c

    γ   = lob.c

    d = ones(T, s)
    ω = get_lobatto_ω_matrix(s)
    δ = zeros(T, 0, s)

    return TableauVSPARKprimary(:vspark_lobabccd, 2*s-2,
                        a_q, a_p, α_q, α_p,
                        a_q̃, a_p̃, α_q̃, α_p̃,
                        b_q, b_p, β_q̃, β_p̃,
                        c_q, c_p, γ, d,
                        ω, δ)
end
