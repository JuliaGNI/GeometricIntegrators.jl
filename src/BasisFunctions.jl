__precompile__()

module BasisFunctions

    export Basis, PolynomialBasis, eval_basis, deriv_basis, int_basis

    include("basis_functions/basis_functions.jl")

    export vandermonde_matrix, vandermonde_matrix_inverse

    include("basis_functions/vandermonde_matrix.jl")

    export BernsteinBasis, BernsteinBasisModified, LagrangeBasis,
           LegendreBasis, LegendreBasisHierarchical

    include("basis_functions/bernstein_basis.jl")
    include("basis_functions/bernstein_basis_modified.jl")
    include("basis_functions/lagrange_basis.jl")
    include("basis_functions/legendre_basis.jl")
    include("basis_functions/legendre_basis_hierarchical.jl")

    export Polynomial, BernsteinPolynomial, LagrangePolynomial

    include("basis_functions/polynomial.jl")


    import ..CommonFunctions

    CommonFunctions.nbasis(b::Basis) = nbasis(b)
    CommonFunctions.nnodes(b::Basis) = nnodes(b)
    CommonFunctions.degree(b::Basis) = degree(b)
    CommonFunctions.nodes(b::Basis)  = nodes(b)

    CommonFunctions.evaluate(b::Basis, x, y) = eval_basis(b, x, y)
    CommonFunctions.derivative(b::Basis, x, y) = deriv_basis(b, x, y)
    CommonFunctions.integral(b::Basis, x, y) = int_basis(b, x, y)

    CommonFunctions.evaluate!(pol::Polynomial, x, y) = evaluate!(pol, x, y)
    CommonFunctions.evaluate!(b::LagrangeBasis, c, x, y) = evaluate!(b, c, x, y)

end
