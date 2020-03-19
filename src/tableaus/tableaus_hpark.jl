
function getTableauHPARK(name, q::CoefficientsRK{T}, p::CoefficientsRK{T}, d=[]) where {T}

    @assert q.s == p.s

    o = min(q.o, p.o)

    a_q = q.a
    a_p = p.a
    b_q = q.b
    b_p = p.b
    c_q = q.c
    c_p = p.c


    if length(d) == 0
        return TableauHPARK(name, o,
                            a_q, a_p, a_q, a_p,
                            a_q, a_p, a_q, a_p,
                            b_q, b_p, b_q, b_p,
                            c_q, c_p, c_q, c_p)
    else
        @assert length(d) == q.s == p.s

        return TableauHPARK(name, o,
                            a_q, a_p, a_q, a_p,
                            a_q, a_p, a_q, a_p,
                            b_q, b_p, b_q, b_p,
                            c_q, c_p, c_q, c_p,
                            d)
    end
end


"Tableau for Gauss-Lobatto IIIA-IIIB HPARK method with two stages."
function getTableauHPARKLobIIIAIIIB2()
    getTableauHPARK(:HPARKLobIIIAIIIB2, getCoefficientsLobIIIA2(), getCoefficientsLobIIIB2())
end

"Tableau for Gauss-Lobatto IIIA-IIIB HPARK method with three stages."
function getTableauHPARKLobIIIAIIIB3()
    getTableauHPARK(:HPARKLobIIIAIIIB3, getCoefficientsLobIIIA3(), getCoefficientsLobIIIB3())
end

"Tableau for Gauss-Lobatto IIIA-IIIB HPARK method with four stages."
function getTableauHPARKLobIIIAIIIB4()
    getTableauHPARK(:HPARKLobIIIAIIIB4, getCoefficientsLobIIIA4(), getCoefficientsLobIIIB4())
end

"Tableau for Gauss-Legendre HPARK method with s stages."
function getTableauHPARKGLRK(s)
    glrk = getCoefficientsGLRK(s)
    getTableauHPARK(Symbol("HPARKGLRK", s), glrk, glrk)
end
