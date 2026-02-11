
abstract type RKMethod <: ODEMethod end
abstract type PRKMethod <: PODEMethod end

const RungeKuttaMethod = Union{RKMethod,PRKMethod}

GeometricBase.tableau(method::RKMethod, args...; kwargs...) = error("No tableau for Runge-Kutta method $(typeof(method)) provided")
GeometricBase.tableau(method::PRKMethod, args...; kwargs...) = error("No tableau for partitioned Runge-Kutta method $(typeof(method)) provided")

GeometricBase.order(method::RungeKuttaMethod) = order(tableau(method))

nstages(method::RungeKuttaMethod) = RungeKutta.nstages(tableau(method))
eachstage(method::RungeKuttaMethod) = RungeKutta.eachstage(tableau(method))
coefficients(method::RungeKuttaMethod) = RungeKutta.coefficients(tableau(method))
weights(method::RungeKuttaMethod) = RungeKutta.weights(tableau(method))
nodes(method::RungeKuttaMethod) = RungeKutta.nodes(tableau(method))

print_reference(io, method::RungeKuttaMethod) =
    try
        ismissing(reference(tableau(method))) || print(io, reference(tableau(method)))
    catch MethodError
        String("")
    end

ispodemethod(::Union{RKMethod,Type{<:RKMethod}}) = true
ishodemethod(::Union{RKMethod,Type{<:RKMethod}}) = true
isiodemethod(::Union{RKMethod,Type{<:RKMethod}}) = true
islodemethod(::Union{RKMethod,Type{<:RKMethod}}) = true
ishodemethod(::Union{PRKMethod,Type{<:PRKMethod}}) = true
isiodemethod(::Union{PRKMethod,Type{<:PRKMethod}}) = true
islodemethod(::Union{PRKMethod,Type{<:PRKMethod}}) = true

isexplicit(method::RungeKuttaMethod) = RungeKutta.isexplicit(tableau(method))
isimplicit(method::RungeKuttaMethod) = RungeKutta.isimplicit(tableau(method))
issymmetric(method::RungeKuttaMethod) = RungeKutta.issymmetric(tableau(method))
issymplectic(method::RungeKuttaMethod) = RungeKutta.issymplectic(tableau(method))
# isenergypreserving(method::RungeKuttaMethod) = RungeKutta.order(tableau(method))
# isstifflyaccurate(method::RungeKuttaMethod) = RungeKutta.order(tableau(method))


function Base.show(io::IO, method::RKMethod)
    print(io, "\nRunge-Kutta Method with Tableau: $(description(tableau(method)))\n")
    print(io, string(tableau(method)))
    print_reference(io, method)
end

function Base.show(io::IO, method::PRKMethod)
    print(io, "\nPartitioned Runge-Kutta Method with Tableau: $(description(tableau(method)))\n")
    print(io, string(tableau(method).q))
    print(io, string(tableau(method).p))
    print_reference(io, method)
end



const StageVector{T} = Vector{<:AbstractVector{T}}

GeometricBase.tableau(int::GeometricIntegrator{<:EquationProblem{ST,DT,TT},<:RKMethod}, args...; kwargs...) where {ST,DT,TT} = tableau(method(int), TT)

initmethod(method::RKMethod, ::GeometricProblem{ST,DT,TT}) where {ST,DT,TT} = RK(method, TT)
initmethod(method::PRKMethod, ::GeometricProblem{ST,DT,TT}) where {ST,DT,TT} = PRK(method, TT)

GeometricIntegrator(problem::AbstractProblemPODE{DT,TT}, method::RKMethod, args...; kwargs...) where {DT,TT} = GeometricIntegrator(problem, PRK(method, TT), args...; kwargs...)
