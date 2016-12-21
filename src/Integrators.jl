__precompile__()

module Integrators

    using ..Equations
    using ..Solutions
    using ..Solvers
    using ..Utils

    import ..Solvers.function_stages!

    include("utils/macro_utils.jl")

    export Tableau, TableauRK, TableauERK, TableauDIRK, TableauFIRK, TableauSIRK,
           TableauEPRK, TableauIPRK, TableauVPRK,
           TableauIPARK, TableauSARK, TableauSPARK,
           TableauGLM,
           showTableau, writeTableauToFile, readTableauERKFromFile

    include("integrators/tableaus.jl")

    export InitialGuess, evaluate, initialize!, update!

    export Integrator, IntegratorERK, IntegratorDIRK, IntegratorFIRK, IntegratorSIRK,
           IntegratorEPRK, IntegratorIPRK, IntegratorVPRK,
           IntegratorIPARK,
           IntegratorSARK, IntegratorSPARK,
           integrate, integrate!, function_stages!

    include("integrators/integrators.jl")
    include("integrators/integrators_erk.jl")
    include("integrators/integrators_eprk.jl")

    include("integrators/initial_guess.jl")

    include("integrators/integrators_dirk.jl")
    include("integrators/integrators_firk.jl")
    include("integrators/integrators_sirk.jl")
    include("integrators/integrators_iprk.jl")
    include("integrators/integrators_vprk.jl")
    include("integrators/integrators_ipark.jl")
    include("integrators/integrators_sark.jl")
    include("integrators/integrators_spark.jl")

end
