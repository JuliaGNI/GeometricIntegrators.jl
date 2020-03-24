module Stochastic

    using LinearAlgebra: dot, mul!
    using OffsetArrays

    using ..CommonFunctions
    using ..Config
    using ..Solvers
    using ..Utils

    import ..Integrators

    import ..Equations: SDE, PSDE, SPSDE, get_function_tuple
    import ..Solutions: AtomicSolutionSDE, AtomicSolutionPSDE, SolutionVector
    import ..Solutions: update!

    import ..Integrators: StochasticIntegrator, Parameters
    import ..Integrators: SDEIntegratorCache, PSDEIntegratorCache,
                          IntegratorCache, CacheDict, CacheType
    import ..Integrators: CoefficientsRK, AbstractTableauERK, AbstractTableauIRK
    import ..Integrators: create_internal_stage_vector, create_internal_stage_matrix,
                          create_internal_stage_vector_with_zero, create_nonlinear_solver


    export StochasticIntegrator, StochasticIntegratorRK, StochasticIntegratorPRK

    export IntegratorSERK, TableauSERK
    export IntegratorSIRK, TableauSIRK
    export IntegratorSIPRK, TableauSIPRK
    export IntegratorSISPRK, TableauSISPRK
    export IntegratorWERK, TableauWERK
    export IntegratorWIRK, TableauWIRK

    include("stochastic/abstract_stochastic_rk.jl")
    include("stochastic/integrators_serk.jl")
    include("stochastic/integrators_sirk.jl")
    include("stochastic/integrators_siprk.jl")
    include("stochastic/integrators_sisprk.jl")
    include("stochastic/integrators_werk.jl")
    include("stochastic/integrators_wirk.jl")
    include("stochastic/common.jl")

end
