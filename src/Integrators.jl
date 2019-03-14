__precompile__()

module Integrators

    using DelimitedFiles
    using Documenter
    using LinearAlgebra
    using OffsetArrays
    using ProgressMeter

    using Base: TwicePrecision

    using ..CommonFunctions
    using ..Config
    using ..Interpolation
    using ..Utils

    using ..BasisFunctions
    using ..Quadratures
    using ..Discontinuities
    using ..Equations
    using ..Solutions
    using ..Solvers


    export InitialGuessODE, InitialGuessPODE, initialize!, update!

    include("integrators/initial_guess/extrapolation.jl")
    include("integrators/initial_guess/initial_guess_ode.jl")
    include("integrators/initial_guess/initial_guess_pode.jl")


    export AbstractCoefficients, AbstractTableau

    include("integrators/abstract_coefficients.jl")
    include("integrators/abstract_tableau.jl")


    export Integrator, integrate, integrate!, equation, timestep
    export NonlinearFunctionParameters, function_stages!

    include("integrators/abstract_integrator.jl")


    include("integrators/integrators_common.jl")


    export CoefficientsRK, AbstractTableauRK, AbstractTableauIRK, AbstractTableauPRK,
           IntegratorRK, writeTableauToFile

    export get_symplectic_conjugate_coefficients, symplecticize,
           check_symplecticity, check_symmetry, compute_symplecticity_error,
           check_order_conditions_B, check_order_conditions_C, check_order_conditions_D

    include("integrators/rk/abstract_integrator_rk.jl")
    include("integrators/rk/coefficients.jl")
    include("integrators/rk/tableaus.jl")


    export IntegratorERK, TableauERK, readTableauERKFromFile

    include("integrators/rk/integrators_erk.jl")


    export IntegratorDIRK, TableauDIRK
    export IntegratorFIRK, TableauFIRK

    include("integrators/rk/integrators_dirk.jl")
    include("integrators/rk/integrators_firk.jl")


    export IntegratorSERK, TableauSERK
    export IntegratorSIRK, TableauSIRK
    export IntegratorSIPRK, TableauSIPRK
    export IntegratorSISPRK, TableauSISPRK
    export IntegratorWERK, TableauWERK
    export IntegratorWIRK, TableauWIRK

    include("integrators/rk/integrators_serk.jl")
    include("integrators/rk/integrators_sirk.jl")
    include("integrators/rk/integrators_siprk.jl")
    include("integrators/rk/integrators_sisprk.jl")
    include("integrators/rk/integrators_werk.jl")
    include("integrators/rk/integrators_wirk.jl")


    export IntegratorEPRK, TableauEPRK
    export IntegratorIPRK, TableauIPRK
    export IntegratorFLRK

    include("integrators/rk/integrators_eprk.jl")
    include("integrators/rk/integrators_iprk.jl")
    include("integrators/rk/integrators_flrk.jl")


    export IntegratorPGLRK, CoefficientsPGLRK

    include("integrators/rk/integrators_pglrk.jl")


    export IntegratorVPRK, IntegratorVPRKpNone, TableauVPRK

    export IntegratorVPRKdegenerate

    export IntegratorVPRKpStandard, IntegratorVPRKpSymplectic,
           IntegratorVPRKpInternal, IntegratorVPRKpMidpoint,
           IntegratorVPRKpSymmetric,
           IntegratorVPRKpSecondary, IntegratorVPRKpVariational

    export IntegratorVPRKpLegendre, TableauVPRKpLegendre

    include("integrators/vprk/integrators_vprk_common.jl")
    include("integrators/vprk/integrators_vprk.jl")
    include("integrators/vprk/integrators_vprk_degenerate.jl")
    include("integrators/vprk/integrators_vprk_pinternal.jl")
    include("integrators/vprk/integrators_vprk_pmidpoint.jl")
    include("integrators/vprk/integrators_vprk_pstandard.jl")
    include("integrators/vprk/integrators_vprk_psecondary.jl")
    include("integrators/vprk/integrators_vprk_psymmetric.jl")
    include("integrators/vprk/integrators_vprk_pvariational.jl")
    include("integrators/vprk/integrators_vprk_plegendre.jl")


    export CoefficientsARK, CoefficientsPRK, CoefficientsMRK

    export TableauARK
    export TableauSARK
    export IntegratorPARK, TableauPARK
    export IntegratorGPARK, TableauGPARK
    export IntegratorVPARK, TableauVPARK
    export IntegratorSPARK, TableauSPARK
    export IntegratorHSPARK, TableauHSPARK
    export IntegratorVSPARK, TableauVSPARK

    include("integrators/spark/coefficients.jl")
    include("integrators/spark/integrators_spark_common.jl")
    include("integrators/spark/integrators_ark.jl")
    include("integrators/spark/integrators_sark.jl")
    include("integrators/spark/integrators_park.jl")
    include("integrators/spark/integrators_gpark.jl")
    include("integrators/spark/integrators_vpark.jl")
    include("integrators/spark/integrators_spark.jl")
    include("integrators/spark/integrators_hspark.jl")
    include("integrators/spark/integrators_vspark.jl")


    export TableauGLM

    include("integrators/glm/integrators_glm.jl")


    export IntegratorSplitting, AbstractTableauSplitting,
           TableauSplittingGS, TableauSplittingNS, TableauSplittingSS

    include("integrators/splitting/integrators_splitting.jl")


    export IntegratorCGVI, IntegratorDGVI, IntegratorDGVIEXP,
           IntegratorDGVIPI, IntegratorDGVIP0, IntegratorDGVIP1

    include("integrators/cgvi/integrators_cgvi.jl")
    include("integrators/dgvi/integrators_dgvi.jl")
    include("integrators/dgvi/integrators_dgvi_experimental.jl")
    include("integrators/dgvi/integrators_dgvi_path_integral.jl")
    include("integrators/dgvi/integrators_dgvi_projection_initial.jl")
    include("integrators/dgvi/integrators_dgvi_projection_final.jl")


    include("integrators/integrators.jl")


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
