module Common

    export State, StateVector

    const State{DT <: Number} = AbstractArray{DT}
    const StateVector{DT,VT} = VT where {DT, VT <: AbstractVector{<:State{DT}}}

    Base.zero(X::ST) where {DT, VT, ST <: StateVector{DT,VT}} = VT[zero(x) for x in X]


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

    export nsamples, ntime, eachsample, eachtimestep

    function nsamples end
    function ntime end
    function eachsample end
    function eachtimestep end

    export nconstraints
    nconstraints() = nothing

end
