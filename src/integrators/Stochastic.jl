module Stochastic

    using LinearAlgebra: dot

    using ..CommonFunctions
    using ..Config
    using ..Solvers
    using ..Utils

    import ..Integrators

    import ..Equations: SDE, PSDE, SPSDE
    import ..Solutions: AtomicSolutionSDE, AtomicSolutionPSDE, SolutionVector

    import ..Integrators: StochasticIntegrator, Parameters
    import ..Integrators: CoefficientsRK, AbstractTableauERK, AbstractTableauIRK
    import ..Integrators: create_internal_stage_vector, create_nonlinear_solver


    export IntegratorSERK, TableauSERK
    export IntegratorSIRK, TableauSIRK
    export IntegratorSIPRK, TableauSIPRK
    export IntegratorSISPRK, TableauSISPRK
    export IntegratorWERK, TableauWERK
    export IntegratorWIRK, TableauWIRK

    include("stochastic/common.jl")
    include("stochastic/integrators_serk.jl")
    include("stochastic/integrators_sirk.jl")
    include("stochastic/integrators_siprk.jl")
    include("stochastic/integrators_sisprk.jl")
    include("stochastic/integrators_werk.jl")
    include("stochastic/integrators_wirk.jl")

end
