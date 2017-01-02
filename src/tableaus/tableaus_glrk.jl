
using FastGaussQuadrature
using Polynomials

function getTableauGLRK(s::Int)

    function evaluate!{T}(pol::Poly{T}, x::Vector{T}, y::Vector{T})
        @assert length(x) == length(y)

        for j in 1:length(y)
            y[j] = pol(x[j])
        end
    end

    # order
    o = 2s

    # obtain Gauss-Legendre nodes and weights
    gl = gausslegendre(s)

    # scale from [-1,+1] to [0,1]
    c = (gl[1]+1)/2
    b = gl[2]/2

    # create Lagrange polynomial
    lag = LagrangePolynomial(c, ones(s))
    vdm = vandermonde_matrix_inverse(c)

    # compute monomial basis functions and corresponding integrals
    polys = []
    poly_ints = []

    for i in 1:s
        y = zeros(s)
        y[i] = 1
        mon = *(vdm, y)
        push!(polys, Poly(mon))
        push!(poly_ints, polyint(polys[i]))
    end

    # compute Runge-Kutta coefficients
    a = zeros(s,s)

    for i in 1:s
        for j in 1:s
            a[i,j] = poly_ints[j](c[i])
        end
    end

    # create tableau
    TableauFIRK(Symbol("glrk", s), o, a, b, c)
end
