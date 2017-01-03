__precompile__()

module BasisFunctions

    export vandermonde_matrix, vandermonde_matrix_inverse

    include("basis_functions/vandermonde_matrix.jl")

    export LagrangePolynomial, evaluate!

    include("basis_functions/lagrange_polynomials.jl")

end
