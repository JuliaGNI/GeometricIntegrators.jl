@doc raw"""
Composition integrator for the solution of initial value problems
```math
\dot{q} (t) = v(t, q(t)) , \qquad q(t_{0}) = q_{0} ,
```
whose vector field ``v`` is given as a sum of vector fields
```math
v (t) = v_1 (t) + ... + v_r (t) .
```

`Composition` has one main constructor:
```julia
Composition(methods, splitting)
```
The first argument is tuple of methods that are used to solve each substep,
and the second argument is a splitting method.
A second convenience constructor uses the [`ExactSolution`](@ref) for all steps:
```julia
Composition(splitting) = Composition(ExactSolution(), splitting)
```
This constructs a composition method that is equivalent to a plain [`Splitting`](@ref).
In order to include exact solutions in the composition, the [`ExactSolution`](@ref) method
implements the general integrator interface.
"""
struct Composition{MT,ST<:AbstractSplittingMethod} <: SODEMethod
    methods::MT
    splitting::ST
end

Composition(splitting::AbstractSplittingMethod) = Composition(ExactSolution(), splitting)

method(c::Composition{<:Tuple}, i::Int) = c.methods[i]
method(c::Composition{<:GeometricMethod}, args...) = c.methods

function methods(c::Composition{<:Tuple}, neqs)
    @assert neqs == length(c.methods)
    return c.methods
end

methods(c::Composition{<:GeometricMethod}, neqs) = Tuple(method(c) for _ in 1:neqs)

splitting(c::Composition) = c.splitting

_options(methods::Tuple) = Tuple([default_options() for _ in methods])
_solvers(methods::Tuple) = Tuple([default_solver(m) for m in methods])
_iguesses(methods::Tuple) = Tuple([default_iguess(m) for m in methods])

_neqs(::Any) = 0
_neqs(::Nothing) = 0
_neqs(eqs::Tuple) = length(eqs)
_neqs(equ::SODE) = nsteps(equ)
_neqs(problem::SODEProblem) = nsteps(problem)


struct CompositionIntegrator{
    MT<:AbstractSplittingMethod,
    PT<:SODEProblem,
    SIT<:Tuple
} <: AbstractIntegrator

    problem::PT
    method::MT
    subints::SIT

    function CompositionIntegrator(
        problem::SODEProblem,
        splitting::AbstractSplittingMethod,
        methods::Tuple;
        options=_options(methods),
        solvers=_solvers(methods),
        initialguesses=_iguesses(methods))

        @assert length(methods) == length(options) == length(solvers) == length(initialguesses) == _neqs(problem)

        # get splitting indices and coefficients
        f, c = coefficients(problem, splitting)

        # construct composition integrators
        subints = Tuple(GeometricIntegrator(SubstepProblem(problem, c[i], f[i]), methods[f[i]]; solver=solvers[f[i]], options[f[i]]...) for i in eachindex(f, c))

        new{typeof(splitting),typeof(problem),typeof(subints)}(problem, splitting, subints)
    end
end

function GeometricIntegrator(problem::SODEProblem, comp::Composition; kwargs...)
    CompositionIntegrator(problem, splitting(comp), methods(comp, _neqs(problem)); kwargs...)
end

problem(int::CompositionIntegrator) = int.problem
subints(int::CompositionIntegrator) = int.subints
method(int::CompositionIntegrator) = int.method

Base.ndims(int::CompositionIntegrator) = ndims(problem(int))

timestep(int::CompositionIntegrator) = timestep(problem(int))

initial_guess!(sol, history, params, ::CompositionIntegrator) = nothing

function initialize!(cint::CompositionIntegrator)
    for int in subints(cint)
        initialize!(int)
    end
end


function integrate_step!(sol, history, params, int::CompositionIntegrator{<:AbstractSplittingMethod,<:SODEProblem})
    # compute composition steps
    for subint in subints(int)
        # compute initial guess for subint
        initial_guess!(sol, history, params, subint)

        # integrate one timestep with subint
        integrate_step!(sol, history, params, subint)
    end
end
