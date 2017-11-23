__precompile__()

module Integrators

    using DoubleDouble
    using ProgressMeter

    using ..CommonFunctions
    using ..Config
    using ..Utils

    using ..BasisFunctions
    using ..Quadratures
    using ..NumericalFluxes
    using ..Equations
    using ..Solutions
    using ..Solvers
    using ..Tableaus


    export InitialGuess, InitialGuessPODE, initialize!, update!

    export Integrator, IntegratorERK, IntegratorDIRK, IntegratorFIRK, IntegratorSIRK,
           IntegratorEPRK, IntegratorIPRK,
           IntegratorPARK, IntegratorVPARK,
           IntegratorSPARK, IntegratorVSPARK,
           IntegratorVPRK, IntegratorVPRKpNone,
           IntegratorVPRKpStandard, IntegratorVPRKpSymplectic,
           IntegratorVPRKpMidpoint, IntegratorVPRKpSymmetric,
           IntegratorVPRKpSecondary, IntegratorVPRKpVariational,
           IntegratorCGVI, IntegratorDGVI,
           IntegratorSplitting,
           integrate, integrate!, function_stages!

    include("integrators/integrators.jl")
    include("integrators/runge_kutta.jl")

    include("integrators/extrapolation.jl")
    include("integrators/initial_guess_ode.jl")
    include("integrators/initial_guess_pode.jl")

    include("integrators/rk/integrators_erk.jl")
    include("integrators/rk/integrators_dirk.jl")
    include("integrators/rk/integrators_firk.jl")
    include("integrators/rk/integrators_sirk.jl")
    include("integrators/rk/integrators_eprk.jl")
    include("integrators/rk/integrators_iprk.jl")

    include("integrators/vprk/integrators_vprk_common.jl")
    include("integrators/vprk/integrators_vprk_pmidpoint.jl")
    include("integrators/vprk/integrators_vprk_pstandard.jl")
    include("integrators/vprk/integrators_vprk_psecondary.jl")
    include("integrators/vprk/integrators_vprk_psymmetric.jl")
    include("integrators/vprk/integrators_vprk_pvariational.jl")
#    include("integrators/vprk/integrators_vprk_pglrk.jl")
    include("integrators/vprk/integrators_vprk.jl")

    include("integrators/spark/integrators_park.jl")
    include("integrators/spark/integrators_vpark.jl")
    include("integrators/spark/integrators_spark.jl")
    include("integrators/spark/integrators_vspark.jl")

    include("integrators/cgvi/integrators_cgvi.jl")
    include("integrators/dgvi/integrators_dgvi.jl")

    include("integrators/integrators_splitting.jl")


    function __init__()
        default_params = (
            (:ig_interpolation, HermiteInterpolation),
            (:ig_extrapolation_stages, 5),
            (:int_show_progress_nmin,  1000),
        )

        for param in default_params
            add_config(param...)
        end
    end

end
