__precompile__()

module BasisFunctions

    export Basis

    include("basis_functions/basis_functions.jl")

    export derivative, integral, evaluate, evaluate!

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

end
