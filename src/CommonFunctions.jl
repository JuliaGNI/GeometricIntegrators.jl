module CommonFunctions

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


    export get_solution, get_solution!, set_solution!, reset!, cut_periodic_solution!

    get_solution() = nothing
    get_solution!() = nothing
    set_solution!() = nothing
    reset!() = nothing
    cut_periodic_solution!() = nothing

end
