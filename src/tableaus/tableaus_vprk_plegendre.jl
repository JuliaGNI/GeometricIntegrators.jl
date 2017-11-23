

"Tableau for variational Gauss-Lobatto IIIA-IIIB method with two stages"
function getTableauVPRKpLegendreQLobIIIA2()
    ω = [0.5, 0.5]
    ω = reshape(ω_λ, 1, length(ω_λ))
    d = [+1.0, -1.0]
    R∞ = -1
    TableauVPRKpLegendreQ(:LobIIIA2, 2, getCoefficientsLobIIIA2(), R∞, ω, d)
end

"Tableau for variational Gauss-Lobatto IIIA-IIIB method with three stages"
function getTableauVPRKpLegendreQLobIIIA3()
    d = [+1.0, -2.0, +1.0]
    R∞ = +1
    TableauVPRKpLegendreQ(:LobIIIA3, 4, getCoefficientsLobIIIA3(), R∞, d)
end

"Tableau for variational Gauss-Lobatto IIIA-IIIB method with four stages"
function getTableauVPRKpLegendreQLobIIIA4()
    d = [+1.0, -√5, +√5, -1.0]
    R∞ = -1
    TableauVPRKpLegendreQ(:LobIIIA4, 6, getCoefficientsLobIIIA4(), R∞, d)
end


"Tableau for variational Gauss-Legendre method with s stages"
function getTableauVPRKpLegendreQGLRK(s; T=Float64)
    glrk = getCoefficientsGLRK(s, T=T)
    R∞ = (-1)^s
    TableauVPRKpLegendreQ(Symbol("GLRK", s), glrk.o, glrk, R∞)
end



"Tableau for variational Gauss-Lobatto IIIA-IIIB method with two stages"
function getTableauVPRKpLegendrePLobIIIA2()
    d = [+1.0, -1.0]
    R∞ = -1
    TableauVPRKpLegendreP(:LobIIIA2, 2, getCoefficientsLobIIIA2(), R∞, d)
end

"Tableau for variational Gauss-Lobatto IIIA-IIIB method with three stages"
function getTableauVPRKpLegendrePLobIIIA3()
    d = [+1.0, -2.0, +1.0]
    R∞ = +1
    TableauVPRKpLegendreP(:LobIIIA3, 4, getCoefficientsLobIIIA3(), R∞, d)
end

"Tableau for variational Gauss-Lobatto IIIA-IIIB method with four stages"
function getTableauVPRKpLegendrePLobIIIA4()
    d = [+1.0, -√5, +√5, -1.0]
    R∞ = -1
    TableauVPRKpLegendreP(:LobIIIA4, 6, getCoefficientsLobIIIA4(), R∞, d)
end


"Tableau for variational Gauss-Legendre method with s stages"
function getTableauVPRKpLegendrePGLRK(s; T=Float64)
    glrk = getCoefficientsGLRK(s, T=T)
    R∞ = (-1)^s
    TableauVPRKpLegendreP(Symbol("GLRK", s), glrk.o, glrk, R∞)
end
