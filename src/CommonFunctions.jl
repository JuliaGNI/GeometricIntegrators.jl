module CommonFunctions

    export State, StateVector

    State{DT} = AbstractArray{DT<:Number}
    StateVector{DT} = AbstractVector{State{DT}}

    export nbasis, nnodes, nodes, order, degree

    nbasis() = nothing
    nnodes() = nothing
    nodes()  = nothing
    order()  = nothing
    degree() = nothing

    export evaluate, evaluate!

    evaluate() = nothing
    evaluate!() = nothing

    export derivative, integral

    derivative() = nothing
    integral() = nothing

    export write_to_hdf5

    write_to_hdf5() = nothing

    export periodicity, reset!, cut_periodic_solution!

    periodicity() = nothing
    reset!() = nothing
    cut_periodic_solution!() = nothing

    export eachsample, eachtimestep

    eachsample() = nothing
    eachtimestep() = nothing

    export nconstraints
    nconstraints() = nothing

end
