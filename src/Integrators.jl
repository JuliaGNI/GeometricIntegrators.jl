__precompile__()

module Integrators

    using DoubleDouble
    using ProgressMeter

    using ..Config
    using ..Utils

    using ..BasisFunctions
    using ..Quadratures
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
    include("integrators/runge_kutta.jl")

    include("integrators/extrapolation.jl")
    include("integrators/initial_guess_ode.jl")
    include("integrators/initial_guess_iode.jl")

    include("integrators/rk/integrators_erk.jl")
    include("integrators/rk/integrators_dirk.jl")
    include("integrators/rk/integrators_firk.jl")
    include("integrators/rk/integrators_sirk.jl")
    include("integrators/rk/integrators_eprk.jl")
    include("integrators/rk/integrators_iprk.jl")

    include("integrators/vprk/integrators_vprk_general.jl")
    include("integrators/vprk/integrators_vprk_pmidpoint.jl")
    include("integrators/vprk/integrators_vprk_pstandard.jl")
    include("integrators/vprk/integrators_vprk_psymplectic.jl")
    include("integrators/vprk/integrators_vprk_psymmetric.jl")
    include("integrators/vprk/integrators_vprk_pvariational.jl")
    include("integrators/vprk/integrators_vprk_pglrk.jl")
    include("integrators/vprk/integrators_vprk.jl")

    include("integrators/spark/integrators_ark.jl")
    include("integrators/spark/integrators_sark.jl")
    include("integrators/spark/integrators_park.jl")
    include("integrators/spark/integrators_spark.jl")
    include("integrators/spark/integrators_vpark.jl")
    include("integrators/spark/integrators_vspark.jl")


    function __init__()
        add_config(:ig_extrapolation_stages, 5)
        add_config(:int_show_progress_nmin,  1000)
    end

end
