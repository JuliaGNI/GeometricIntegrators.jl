@doc raw"""
Splitting integrator for the solution of initial value problems
```math
\dot{q} (t) = v(t, q(t)) , \qquad q(t_{0}) = q_{0} ,
```
whose vector field ``v`` is given as a sum of vector fields
```math
v (t) = v_1 (t) + ... + v_r (t) .
```

`Splitting` has two constructors:
```julia
IntegratorSplitting{DT,D}(solutions::Tuple, f::Vector{Int}, c::Vector, Δt)
IntegratorSplitting(equation::SODE, tableau::AbstractTableauSplitting, Δt)
```
In the first constructor, `DT` is the data type of the state vector and `D`
the dimension of the system. In the second constructor, this information
is extracted from the equation. 
The tuple `solutions` contains functions implementing the flow (exact solution)
of the vector fields `v_i`. The vectors `f` and `c` define the actual splitting
method: `f` is a vector of indices of the flows in the split equation to be
solved and `c` is a vector of the same size `f` that contains the coefficients
for each splitting step, i.e., the resulting integrator has the form
```math
\varphi_{\tau} = \phi_{c[s] \tau}^{v_{f[s]}} \circ \dotsc \circ \phi_{c[2] \tau}^{v_{f[2]}} \circ \phi_{c[1] \tau}^{v_{f[1]}} .
```
In the second constructor, these vectors are constructed from the tableau and
the equation.

"""
struct Splitting{T} <: SODEMethod
    f::Vector{Int}
    c::Vector{T}
end

Splitting(method::AbstractSplittingMethod, problem::SODEProblem) = Splitting(coefficients(problem, method)...)

coefficients(method::Splitting) = (method.f, method.c)

initmethod(method::AbstractSplittingMethod, problem::SODEProblem) = Splitting(method, problem)
initmethod(method::Splitting, ::SODEProblem) = method


"Splitting integrator cache."
mutable struct IntegratorCacheSplitting{DT,TT,D,AT} <: ODEIntegratorCache{DT,D}
    q::AT
    t::TT

    function IntegratorCacheSplitting{DT,TT,D}(q₀::AT) where {DT,TT,D,AT <: AbstractArray{DT}}
        new{DT,TT,D,AT}(zero(q₀), zero(TT))
    end
end

function Cache{ST}(problem::SODEProblem, method::Splitting; kwargs...) where {ST}
    IntegratorCacheSplitting{ST, typeof(timestep(problem)), ndims(problem)}(problem.ics.q; kwargs...)
end

@inline CacheType(ST, problem::SODEProblem, ::Splitting) = IntegratorCacheSplitting{ST, typeof(timestep(problem)), ndims(problem)}


function integrate_step!(
    solstep::SolutionStepODE{DT,TT},
    problem::SODEProblem{DT,TT},
    method::Splitting,
    caches::CacheDict,
    ::NoSolver) where {DT,TT}
    
    # compute splitting steps
    for i in eachindex(method.f, method.c)
        if method.c[i] ≠ zero(TT)
            # copy previous solution
            caches[DT].q .= solstep.q

            # compute new solution
            solutions(problem).q[method.f[i]](solstep.q, solstep.t̄[1] + timestep(problem) * method.c[i], caches[DT].q, solstep.t̄[1])
        end
    end
end
