
using FastGaussQuadrature


function getCoefficientsPGLRK(s::Int)
    @assert s ≥ 2

    # order
    o = 2s

    # obtain Gauss-Legendre nodes and weights
    gl = gausslegendre(s)

    # scale from [-1,+1] to [0,1]
    c = (gl[1] .+ 1) ./ 2
    b = gl[2] ./ 2

    T = eltype(c)

    local ξᵢ::T

    leg_basis = LegendreBasis(T,s)

    a = zeros(T, s, s)
    t = zeros(T, s, s)
    P = zeros(T, s, s)
    X = zeros(T, s, s)
    W = zeros(T, s, s)

    for j in 1:s
        for i in 1:s
            P[i,j] = evaluate(leg_basis, j, c[i])
        end
    end

    Q = inv(P)

    X[1,1] = 1/2

    for i in 1:s-1
        ξᵢ = one(T)/sqrt(convert(T, 4i^2-1))/2
        X[i+1,i] = +ξᵢ
        X[i,i+1] = -ξᵢ
    end

    W[s, s-1] = +1
    W[s-1, s] = -1

    mul!(t, X, Q)
    mul!(a, P, t)

    CoefficientsPGLRK(Symbol("pglrk", s), o, a, b, c, P, X, W)
end
