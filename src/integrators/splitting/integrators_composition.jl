@doc raw"""
Composition integrator for the solution of initial value problems
```math
\dot{q} (t) = v(t, q(t)) , \qquad q(t_{0}) = q_{0} ,
```
whose vector field ``v`` is given as a sum of vector fields
```math
v (t) = v_1 (t) + ... + v_r (t) .
```

`IntegratorComposition` has three constructors:
```julia
IntegratorComposition{DT,D}(integrators::Tuple, Δt)
IntegratorComposition(equation::SODE, constructors::Tuple, tableau::AbstractTableauSplitting, Δt)
IntegratorComposition(equation::SODE, tableau::AbstractTableauSplitting, Δt)
```
In the first constructor, `DT` is the data type of the state vector and `D`
the dimension of the system. In the second and third constructor, this
information is extracted from the equation. 
The tuple `integrators` contains the integrators for each substep. Each integrator
is instantiated with appropriately scaled time step size $\Delta t = c_i \tau$ to
match the corresponding splitting scheme.
In the second constructor, the tuple `constructors` contains constructors for the
integrators of each step of the composition. The integrators are constructed
according to the tableau and time step `\Delta t` and passed to the first
constructor.
The third constructor assumes that the exact solution is used for each splitting
step. It thus constructs a composition method that is equivalent to a plain
[`IntegratorSplitting`](@ref).

In order to include exact solutions in the composition, the [`IntegratorExactODE`](@ref)
implements the general integrator interface.

"""
struct Composition{T,IT}
    f::Vector{Int}
    c::Vector{T}
    integrators::IT
end


function Composition(problem::SODEProblem,  splitting::SplittingCoefficients, integrators::Tuple)
    _functions = functions(problem).v
    _solutions = solutions(problem).q

    D = ndims(problem)
    R = length(_functions)

    f, c = splitting_coefficients(R, splitting)

    subints = ()

    # construct composition integrators
    for i in eachindex(f,c)
        if c[i] ≠ zero(TT)
            cᵢ = Δt * c[i]

            if integrators[f[i]] == IntegratorExactODE
                @assert hassolution(equation(problem), f[i])
                subint = integrators[f[i]]{DT,D}(_solutions[f[i]], cᵢ)
            else
                @assert hasvectorfield(equation(problem), f[i])
                subint = integrators[f[i]](_functions[f[i]], cᵢ)
            end

            subints  = (subints..., subint)
        end
    end

    Composition(f, c, subints)
end


function Composition(problem::SODEProblem{DT}, splitting::SplittingCoefficients) where {DT}
    @assert hassolution(equation(problem))
    solutions = Tuple(ExactSolutionODE(q) for (c,q) in zip(coefficients(problem, splitting), solutions(problem).q))
    Composition(problem, solutions, coeffs, timestep(problem))
end


function Cache{ST}(problem::SODEProblem, method::Composition; kwargs...) where {ST}
    IntegratorCacheSplitting{ST, typeof(timestep(problem)), ndims(problem)}(; kwargs...)
end

@inline CacheType(ST, problem::SODEProblem, ::Composition) = IntegratorCacheSplitting{ST, typeof(timestep(problem)), ndims(problem)}


function integrate_step!(
    solstep::SolutionStepODE{DT,TT},
    problem::SODEProblem{DT,TT},
    method::Composition,
    caches::CacheDict,
    solver::NoSolver) where {DT,TT}
    
    # compute composition steps
    for subint in method.integrators
        # copy previous solution
        caches[DT].t .= solstep.t
        caches[DT].q .= solstep.q

        # initialize!(subint, sol)

        integrate_step!(solstep, problem, subint, caches, solver)
    end
end
