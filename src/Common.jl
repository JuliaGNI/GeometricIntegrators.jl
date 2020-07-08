module Common

    export State, StateVector

    State{DT <: Number} = AbstractArray{DT}
    StateVector{DT} = AbstractVector{<:State{DT}}

    Base.zero(X::AT) where {AT <: StateVector} = AT[zero(x) for x in X]


    export nbasis, nnodes, nodes, order, degree

    function nbasis end
    function nnodes end
    function nodes end
    function order end
    function degree end

    export evaluate, evaluate!

    function evaluate end
    function evaluate! end

    export derivative, integral

    function derivative end
    function integral end

    export write_to_hdf5

    function write_to_hdf5 end

    export periodicity, reset!, cut_periodic_solution!

    function periodicity end
    function reset! end
    function cut_periodic_solution! end

    export eachsample, eachtimestep

    function eachsample end
    function eachtimestep end

end
