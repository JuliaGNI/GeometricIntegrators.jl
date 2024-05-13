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
Splitting(f::Vector{Int}, c::Vector)
Splitting(method::AbstractSplittingMethod, problem::SODEProblem)
```
The vectors `f` and `c` define the actual splitting method: `f` is a vector
of indices of the flows in the split equation (that is the exact solution which
can be obtained by `solutions(problem)`) to be solved
and `c` is a vector of the same size as `f` that contains the coefficients
for each splitting step, i.e., the resulting integrator has the form
```math
\varphi_{\tau} = \phi_{c[s] \tau}^{v_{f[s]}} \circ \dotsc \circ \phi_{c[2] \tau}^{v_{f[2]}} \circ \phi_{c[1] \tau}^{v_{f[1]}} .
```
In the second constructor, these vectors are constructed from the splitting method and
the problem.
"""
struct Splitting{T} <: SODEMethod
    f::Vector{Int}
    c::Vector{T}
end

Splitting(method::AbstractSplittingMethod, problem::SODEProblem) = Splitting(coefficients(problem, method)...)
Splitting(problem::SODEProblem, method::AbstractSplittingMethod) = Splitting(coefficients(problem, method)...)

coefficients(method::Splitting) = (f = method.f, c = method.c)

initmethod(method::AbstractSplittingMethod, problem::SODEProblem) = Splitting(method, problem)
initmethod(method::Splitting, ::SODEProblem) = method


"Splitting integrator cache."
mutable struct SplittingCache{DT,TT,D,AT} <: ODEIntegratorCache{DT,D}
    q::AT
    t::TT

    function SplittingCache{DT,TT,D}(q₀::AT) where {DT,TT,D,AT <: AbstractArray{DT}}
        new{DT,TT,D,AT}(zero(q₀), zero(TT))
    end
end

function Cache{ST}(problem::SODEProblem, method::Splitting; kwargs...) where {ST}
    SplittingCache{ST, typeof(timestep(problem)), ndims(problem)}(initial_conditions(problem).q; kwargs...)
end

@inline CacheType(ST, problem::SODEProblem, ::Splitting) = SplittingCache{ST, typeof(timestep(problem)), ndims(problem)}

function reset!(cache::SplittingCache, t, q, λ = missing)
    copyto!(cache.q, q)
    cache.t = t
end

function integrate_step!(sol, history, params, int::GeometricIntegrator{<:Splitting, <:SODEProblem})
    # compute splitting steps
    for i in eachindex(method(int).f, method(int).c)
        if method(int).c[i] ≠ 0
            # copy previous solution and compute time
            cache(int).q .= sol.q
            cache(int).t  = sol.t + timestep(int) * (method(int).c[i] - 1)

            # compute new solution
            solutions(problem(int)).q[method(int).f[i]](sol.q, cache(int).t, cache(int).q, sol.t - timestep(int), params)
        end
    end
end

function integrate_step!(int::GeometricIntegrator{<:Splitting, <:SODEProblem})
    integrate_step!(current(solstep(int)), history(solstep(int)), parameters(solstep(int)), int)
end
