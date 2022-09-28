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
struct IntegratorComposition{DT, TT, D, S, IT <: Tuple} <: ODEIntegrator{DT,TT}
    ints::IT
    Δt::TT
    q̄::Vector{DT}

    function IntegratorComposition{DT,D}(integrators::IT, Δt::TT) where {DT, TT, D, IT <: Tuple}
        new{DT, TT, D, length(integrators), IT}(integrators, Δt, zeros(DT,D))
    end
end

function IntegratorComposition(problem::SODEProblem{DT,TT}, integrators::Tuple, tableau::ST, Δt::TT=tstep(problem)) where {DT, TT, ST <: AbstractTableauSplitting{TT}}
    _functions = functions(problem).v
    _solutions = solutions(problem).q

    D = ndims(problem)
    R = length(_functions)

    f, c = get_splitting_coefficients(R, tableau)

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

    IntegratorComposition{DT, ndims(problem)}(subints, Δt)
end

function IntegratorComposition(problem::SODEProblem{DT}, tableau::AbstractTableauSplitting) where {DT}
    @assert hassolution(equation(problem))
    # integrators = Tuple(IntegratorConstructor(DT,ndims(equation)) for q in solutions(problem).q)
    integrators = Tuple(IntegratorExactODE for q in solutions(problem).q)
    IntegratorComposition(problem, integrators, tableau, timestep(problem))
end


@inline Base.ndims(::IntegratorComposition{DT,TT,D}) where {DT,TT,D} = D
timestep(int::IntegratorComposition) = int.Δt


function integrate_step!(int::IntegratorComposition{DT,TT}, sol::SolutionStepODE{DT,TT}) where {DT,TT}
    local t̄ = sol.t
    int.q̄  .= sol.q

    # compute composition steps
    for subint in int.ints
        initialize!(subint, sol)
        integrate_step!(subint, sol)
        sol.t̄  = t̄
        sol.q̄ .= int.q̄
    end
end
