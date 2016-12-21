__precompile__()

module Tableaus

    using ..Integrators

    include("tableaus/tableaus_erk.jl")

    export getTableauExplicitEuler, getTableauExplicitMidpoint, getTableauHeun,
           getTableauKutta, getTableauERK4, getTableauERK438

    include("tableaus/tableaus_dirk.jl")

    export getTableauCrouzeix

    include("tableaus/tableaus_firk.jl")

    export getTableauImplicitEuler, getTableauImplicitMidpoint,
           getTableauGLRK1, getTableauGLRK2, getTableauGLRK3

    include("tableaus/tableaus_sirk.jl")

    include("tableaus/tableaus_eprk.jl")

    export getTableauSymplecticEulerA, getTableauSymplecticEulerB

    include("tableaus/tableaus_vprk.jl")

    export getTableauLobIIIAB2

    include("tableaus/tableaus_ipark.jl")

    export getTableauSymplecticProjection

    include("tableaus/tableaus_spark.jl")

end
