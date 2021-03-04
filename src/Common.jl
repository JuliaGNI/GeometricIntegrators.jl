module Common

    export OptionalAbstractArray, OptionalArray, OptionalFunction, OptionalNamedTuple
    
    const OptionalArray{arrayType} = Union{Nothing, arrayType} where {arrayType <: AbstractArray}
    
    const OptionalAbstractArray = Union{Nothing, AbstractArray}
    const OptionalFunction      = Union{Nothing, Function}
    const OptionalNamedTuple    = Union{Nothing, NamedTuple}


    export State, StateVector

    const State{DT <: Number} = AbstractArray{DT}
    const StateVector{DT,VT} = VT where {DT, VT <: AbstractVector{<:State{DT}}}

    Base.zero(X::ST) where {DT, VT, ST <: StateVector{DT,VT}} = VT[zero(x) for x in X]


    export evaluate, evaluate!

    function evaluate end
    function evaluate! end

    export write_to_hdf5

    function write_to_hdf5 end

    export periodicity, reset!, cut_periodic_solution!

    function periodicity end
    function reset! end
    function cut_periodic_solution! end

    export nsamples, nconstraints, ntime, eachsample, eachtimestep

    function nsamples end
    function nconstraints end
    function ntime end
    function eachsample end
    function eachtimestep end

end
