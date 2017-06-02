__precompile__()

module BasisFunctions

    export vandermonde_matrix, vandermonde_matrix_inverse

    include("basis_functions/vandermonde_matrix.jl")

    export LagrangeBasis, LagrangePolynomial, lagrange, lagrange_derivative, evaluate!

    include("basis_functions/lagrange_polynomials.jl")

end
