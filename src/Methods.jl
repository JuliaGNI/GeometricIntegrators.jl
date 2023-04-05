module Methods

    using GeometricBase
    using GeometricBase.Utils
    using GeometricEquations
    using PrettyTables
    using RungeKutta.Tableaus
    using RungeKutta.PartitionedTableaus

    import RungeKutta
    import RungeKutta: AbstractTableau, Tableau, PartitionedTableau, SymplecticTableau, SymplecticPartitionedTableau
    import RungeKutta: eachstage, nstages
    import SimpleSolvers: SolverMethod


    include("methods/methods.jl")
    include("methods/rungekutta.jl")
    # include("methods/splitting.jl")
    include("methods/flrk.jl")
    include("methods/vprk.jl")
    include("methods/vprk_degenerate.jl")
    include("methods/dvi.jl")

    include("methods/projection.jl")
    include("methods/vprk_projected.jl")

    export GeometricMethod
    export ODEMethod, PODEMethod, HODEMethod, IODEMethod, LODEMethod, SODEMethod
    export DAEMethod, PDAEMethod, HDAEMethod, IDAEMethod, LDAEMethod
    export RKMethod, ERKMethod, IRKMethod, DIRKMethod
    export PRKMethod, EPRKMethod, IPRKMethod, VPRKMethod
    export AbstractSplittingMethod

    export initmethod

    export RK, PRK
    export internal_variables
    export implicit_update

    export NoProjection, projection

    export AbstractTableau, Tableau, PartitionedTableau, SymplecticTableau, SymplecticPartitionedTableau
    export nstages, eachstage, coefficients, weights, nodes
    export hasnullvector, nullvector

    include("methods/list.jl")

end
