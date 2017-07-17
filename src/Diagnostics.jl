__precompile__()

module Diagnostics

    using HDF5

    using ..CommonFunctions
    using ..BasisFunctions
    using ..Equations
    using ..Integrators
    using ..Solutions
    using ..Tableaus
    using ..Utils

    export PoincareInvariant1st, PoincareInvariant2nd, PoincareInvariant2ndTrapezoidal,
           evaluate_poincare_invariant, write_to_hdf5

    include("diagnostics/poincare_invariant_1st.jl")
    include("diagnostics/poincare_invariant_2nd.jl")
    include("diagnostics/poincare_invariant_2nd_trapezoidal.jl")

end
