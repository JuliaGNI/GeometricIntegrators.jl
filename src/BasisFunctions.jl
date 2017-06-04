__precompile__()

module BasisFunctions

    export derivative, integral

    export vandermonde_matrix, vandermonde_matrix_inverse

    include("basis_functions/vandermonde_matrix.jl")

    export LagrangeBasis, evaluate

    include("basis_functions/lagrange_basis.jl")

    export LagrangePolynomial, evaluate!

    include("basis_functions/lagrange_polynomials.jl")

end
