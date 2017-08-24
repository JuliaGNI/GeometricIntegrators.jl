__precompile__()

module Diagnostics

    using HDF5
    using ProgressMeter

    using ..CommonFunctions
    using ..Config
    using ..Utils

    using ..BasisFunctions
    using ..Equations
    using ..Integrators
    using ..Solutions
    using ..Tableaus

    export PoincareInvariant1st, PoincareInvariant1stCanonical,
           PoincareInvariant2nd, PoincareInvariant2ndCanonical,
           PoincareInvariant2ndTrapezoidal,
           evaluate_poincare_invariant, write_to_hdf5

    include("diagnostics/poincare_invariant_1st_common.jl")
    include("diagnostics/poincare_invariant_1st.jl")
    include("diagnostics/poincare_invariant_1st_canonical.jl")
    include("diagnostics/poincare_invariant_2nd.jl")
    include("diagnostics/poincare_invariant_2nd_canonical.jl")
    include("diagnostics/poincare_invariant_2nd_trapezoidal.jl")
    include("diagnostics/poincare_invariant_2nd_common.jl")

end
