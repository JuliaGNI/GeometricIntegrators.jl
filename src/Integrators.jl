__precompile__()

module Integrators

    using ..Config
    using ..Utils
    using ..Equations
    using ..Solutions
    using ..Solvers
    using ..Tableaus

    import ..Solvers.function_stages!

    export InitialGuess, InitialGuessIODE, initialize!, update!

    export Integrator, IntegratorERK, IntegratorDIRK, IntegratorFIRK, IntegratorSIRK,
           IntegratorEPRK, IntegratorIPRK, IntegratorVPRK, IntegratorVPARK, IntegratorVSPARK,
           IntegratorVPRKpStandard, IntegratorVPRKpSymplectic, IntegratorVPRKpSymmetric,
           IntegratorVPRKpMidpoint, IntegratorIPARK,
           IntegratorSARK, IntegratorSPARK,
           integrate, integrate!, function_stages!

    include("integrators/integrators.jl")
    include("integrators/integrators_erk.jl")

    include("integrators/initial_guess_ode.jl")
    include("integrators/initial_guess_iode.jl")

    include("integrators/integrators_dirk.jl")
    include("integrators/integrators_firk.jl")
    include("integrators/integrators_sirk.jl")
    include("integrators/integrators_eprk.jl")
    include("integrators/integrators_iprk.jl")
    include("integrators/integrators_vprk_general.jl")
    include("integrators/integrators_vprk_pmidpoint.jl")
    include("integrators/integrators_vprk_pstandard.jl")
    include("integrators/integrators_vprk_psymmetric.jl")
    include("integrators/integrators_vprk.jl")
    include("integrators/integrators_ark.jl")
    include("integrators/integrators_sark.jl")
    include("integrators/integrators_park.jl")
    include("integrators/integrators_spark.jl")
    include("integrators/integrators_vpark.jl")
    include("integrators/integrators_vspark.jl")

end
