
using .SPARK
using .VPRK

#*****************************************************************************#
# General initialization functions for all integrators                        #
#*****************************************************************************#

"""
```julia
Integrator(equation, tableau, Δt)
```

Universal constructor for Runge-Kutta (standard, partitioned, variational),
SPARK and splitting integrators that automatically selects the appropriate
integrator based on the equation and tableau types.

"""
function Integrator end


# Print error for integrators not implemented, yet.
function Integrator(problem::GeometricProblem, tableau::Union{AbstractTableau,Tableau})
    error("No integrator found for problem type ", typeof(problem), " and tableau ", tableau)
end


#*****************************************************************************#
# Initialization functions for deterministic integrators                      #
#*****************************************************************************#

# Create integrator for Runge-Kutta tableau.
function Integrator(problem::ODEProblem, tableau::Tableau)
    if isexplicit(tableau)
        # Create integrator for explicit Runge-Kutta tableau
        IntegratorERK(problem, tableau)
    elseif isdiagnonallyimplicit(tableau)
        # Create integrator for diagonally implicit Runge-Kutta tableau
        IntegratorDIRK(problem, tableau)
    elseif isfullyimplicit(tableau)
        # Create integrator for fully implicit Runge-Kutta tableau
        IntegratorFIRK(problem, tableau)
    end
end

# Create integrator for explicit partitioned Runge-Kutta tableau.
function Integrator(problem::Union{PODEProblem,HODEProblem}, tableau::PartitionedTableau)
    if isexplicit(tableau)
        IntegratorEPRK(problem, tableau)
    else
        IntegratorIPRK(problem, tableau)
    end
end

function Integrator(problem::Union{PODEProblem,HODEProblem}, tableau::Tableau)
    Integrator(problem, PartitionedTableau(tableau))
end

# Create integrator for implicit partitioned Runge-Kutta tableau.
function Integrator(problem::Union{IODEProblem,LODEProblem}, tableau::PartitionedTableau)
    IntegratorPRKimplicit(problem, tableau)
end

# Create integrator for variational partitioned Runge-Kutta tableau.
function Integrator(problem::Union{IODEProblem,LODEProblem}, tableau::TableauVPRK)
    IntegratorVPRK(problem, tableau)
end

# Create integrator for formal Lagrangian Runge-Kutta tableau.
function Integrator(problem::LODEProblem, tableau::Tableau)
    IntegratorFLRK(problem, tableau)
end

# Create integrator for Projected Gauss-Legendre Runge-Kutta tableau.
function Integrator(problem::Union{IODEProblem,LODEProblem}, tableau::CoefficientsPGLRK)
    IntegratorPGLRK(problem, tableau)
end

# Create integrator for variational partitioned additive Runge-Kutta tableau.
function Integrator(problem::Union{IDAEProblem,LDAEProblem}, tableau::TableauVPARK)
    IntegratorVPARK(problem, tableau)
end

# Create integrator for special partitioned additive Runge-Kutta tableau.
function Integrator(problem::Union{IDAEProblem,LDAEProblem}, tableau::TableauSPARK)
    IntegratorSPARK(problem, tableau)
end

# Create integrator for variational special partitioned additive Runge-Kutta tableau.
function Integrator(problem::Union{IDAEProblem,LDAEProblem}, tableau::TableauVSPARK)
    IntegratorVSPARK(problem, tableau)
end

# Create integrator for variational special partitioned additive Runge-Kutta tableau with projection on primary constraint.
function Integrator(problem::Union{IDAEProblem,LDAEProblem}, tableau::TableauVSPARKprimary)
    IntegratorVSPARKprimary(problem, tableau)
end

# Create integrator for variational special partitioned additive Runge-Kutta tableau with projection on secondary constraint.
function Integrator(problem::LDAEProblem, tableau::TableauVSPARKsecondary)
    IntegratorVSPARKsecondary(problem, tableau)
end

# Create integrator for Hamiltonian partitioned additive Runge-Kutta tableau.
function Integrator(problem::Union{PDAEProblem,HDAEProblem}, tableau::TableauHPARK)
    IntegratorHPARK(problem, tableau)
end

# Create integrator for Hamiltonian special partitioned additive Runge-Kutta tableau.
function Integrator(problem::Union{PDAEProblem,HDAEProblem}, tableau::TableauHSPARK)
    IntegratorHSPARK(problem, tableau)
end

# Create integrator for Hamiltonian special partitioned additive Runge-Kutta tableau with projection on primary constraint.
function Integrator(problem::Union{PDAEProblem,HDAEProblem}, tableau::TableauHSPARKprimary)
    IntegratorHSPARKprimary(problem, tableau)
end

# Create integrator for splitting tableau.
function Integrator(problem::SODEProblem, tableau::AbstractTableauSplitting)
    IntegratorSplitting(problem, tableau)
end


#*****************************************************************************#
# Constructor wrappers for deterministic integrators                          #
#*****************************************************************************#

"""
```julia
IntegratorConstructor(DT, D)
IntegratorConstructor(DT, D, tableau)
```

Returns a closure for various integrator constructors that is needed for
composition methods based on the tableau type. If not tableau is provided,
a closure for the constructor of an exact solution is returned.

- `DT`: data type of the solution / state vector
- `D`: dimension of the state vector

"""
function IntegratorConstructor end

# Create integrator constructor for exact solution.
function IntegratorConstructor(DT, D)
    (v::Function, Δt::Number; kwargs...) -> IntegratorExactODE{DT,D}(v, Δt; kwargs...)
end

# Create integrator constructor for Runge-Kutta tableau.
function IntegratorConstructor(DT, D, tableau::Tableau)
    if isexplicit(tableau)
        # Create integrator constructor for explicit Runge-Kutta tableau
        (v::Function, Δt::Number; kwargs...) -> IntegratorERK{DT,D}(v, tableau, Δt; kwargs...)
    elseif isdiagnonallyimplicit(tableau)
        # Create integrator constructor for diagonally implicit Runge-Kutta tableau
        (v::Function, Δt::Number; kwargs...) -> IntegratorDIRK{DT,D}(v, tableau, Δt; kwargs...)
    elseif isfullyimplicit(tableau)
        # Create integrator constructor for fully implicit Runge-Kutta tableau
        (v::Function, Δt::Number; kwargs...) -> IntegratorFIRK{DT,D}(v, tableau, Δt; kwargs...)
    end
end


#*****************************************************************************#
# General integration functions for all integrators                           #
#*****************************************************************************#

"""
```julia
integrate(equation, integrator, ntime; kwargs...)
integrate(equation, tableau, Δt, ntime; kwargs...)
```

Integrate an `equation` with `integrator` for `ntime` time steps
and return the solution.
If a `tableau` and a time step `Δt` is passed instead of an
`integrator`, the appropriate integrator is created automatically.

Some convenience methods exist for the integration of ODEs,
```julia
integrate(v, x₀, tableau, Δt, ntime; t₀=0., kwargs...)
```
where `v` is the function for the vector field, `x₀` the initial condition
and `t₀` the initial time, and for PODEs
```julia
integrate(v, f, q₀, p₀, tableau, Δt, ntime; t₀=0., kwargs...)
```
with vector fields `v` and `f` and initial conditions `q₀` and `p₀`.
"""
function integrate end

# Apply integrator for ntime time steps and return solution.
function integrate(problem::GeometricProblem, integrator::Integrator; kwargs...)
    solution = Solution(problem; kwargs...)
    integrate!(integrator, solution)
    return solution
end

# Integrate given equation with given tableau for ntime time steps and return solution.
function integrate(problem::GeometricProblem, tableau::Union{AbstractTableau,Tableau,PartitionedTableau}; kwargs...)
    return integrate(problem, Integrator(problem, tableau); kwargs...)
end

# Integrate ODE specified by vector field and initial condition with given tableau for ntime time steps and return solution.
function integrate(f::Function, x₀::Vector, tableau::Union{AbstractTableau,Tableau}, tspan, Δt; kwargs...)
    return integrate(ODEProblem(f, tspan, Δt, x₀), tableau; kwargs...)
end

# Integrate PODE specified by two vector fields and initial conditions with given tableau for ntime time steps and return solution.
function integrate(v::Function, f::Function, q₀::Vector, p₀::Vector, tableau::Union{AbstractTableau,Tableau,PartitionedTableau}, tspan, Δt; kwargs...)
    return integrate(PODEProblem(v, f, tspan, Δt, q₀, p₀), tableau; kwargs...)
end


#*****************************************************************************#
# Integration functions for deterministic integrators                         #
#*****************************************************************************#

"""
```julia
integrate!(integrator, solution)
```

Perform the actual integration for a given `integrator` and `solution`.
Solve all time steps for all initial conditions.

```julia
integrate!(integrator, solution, m₁, m₂)
```

Solve all time steps for initial conditions m with m₁ ≤ m ≤ m₂.

```julia
integrate!(integrator, solution, m₁, m₂, n₁, n₂)
```

Solve time steps n with n₁ ≤ n ≤ n₂ for initial conditions m with m₁ ≤ m ≤ m₂.

```julia
integrate!(integrator, solution, atomic_solution, m, n)
```

Solve one time step n for one initial condition m.

"""
function integrate! end

# Parts of one integration step that are common to deterministic and stochastic equations.
function integrate_common!(int::Integrator{DT,TT}, sol::AbstractSolution{AT,TT}, asol::AtomicSolution{DT,TT}, m::Int, n::Int) where {DT, TT, AT <: AbstractArray{DT}}
    # integrate one initial condition for one time step
    integrate_step!(int, asol)

    # take care of periodic solutions
    cut_periodic_solution!(asol, periodicity(sol))

    # copy solution from cache to solution
    set_solution!(sol, asol, n, m)
end

# Integrate one time step n for one initial condition m.
function integrate!(int::Integrator{DT,TT}, sol::AbstractSolution{AT,TT}, asol::AtomicSolution{DT,TT}, m::Int, n::Int) where {DT, TT, AT <: AbstractArray{DT}}
    integrate_common!(int, sol, asol, m, n)
end

# Integrate `equation` for all initial conditions and all time steps.
function integrate!(int::Integrator, sol::AbstractSolution)
    integrate!(int, sol, 1, nsamples(sol))
end


# Integrate equation for initial conditions m with m₁ ≤ m ≤ m₂.
function integrate!(int::Integrator, sol::AbstractSolution, m₁, m₂)
    # integrate samples m with m₁ ≤ m ≤ m₂ for all time steps
    integrate!(int, sol, m₁, m₂, 1, ntime(sol))
end


# Integrate equation for initial conditions m with m₁ ≤ m ≤ m₂ for time steps n with n₁ ≤ n ≤ n₂.
function integrate!(int::Integrator{DT,TT}, sol::AbstractSolution{AT,TT}, m₁::Int, m₂::Int, n₁::Int, n₂::Int) where {DT, TT, AT <: AbstractArray{DT}}
    @assert m₁ ≥ 1
    @assert m₂ ≥ m₁
    @assert m₂ ≤ nsamples(sol)

    # TODO: reactivate after fixing TimeSeries
    # @assert n₁ ≥ 1
    # @assert n₂ ≥ n₁
    # @assert n₂ ≤ ntime(sol)

    asol = AtomicSolution(sol, int)

    # loop over initial conditions showing progress bar
    for m in m₁:m₂
        # get cache from solution
        get_initial_conditions!(sol, asol, m, n₁)
        initialize!(int, asol)

        # loop over time steps
        for n in n₁:n₂
            integrate!(int, sol, asol, m, n)

            # try
            #     integrate!(int, sol, asol, m, n)
            # catch ex
            #     tstr = " in time step " * string(n)
            #
            #     if m₁ ≠ m₂
            #         tstr *= " for initial condition " * string(m)
            #     end
            #
            #     tstr *= "."
            #
            #     if isa(ex, DomainError)
            #         @warn("Domain error" * tstr)
            #     elseif isa(ex, ErrorException)
            #         @warn("Simulation exited early" * tstr)
            #         @warn(ex.msg)
            #     else
            #         @warn(string(typeof(ex)) * tstr)
            #         throw(ex)
            #     end
            # end
        end
    end
end
