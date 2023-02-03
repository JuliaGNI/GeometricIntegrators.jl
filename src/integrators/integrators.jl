
using .SPARK
using .VPRK

#*****************************************************************************#
# General initialization functions for all integrators                        #
#*****************************************************************************#

"""
```julia
Integrator(problem, tableau)
```

Universal constructor for Runge-Kutta (standard, partitioned, variational),
SPARK and splitting integrators that automatically selects the appropriate
integrator based on the problem and tableau types.
"""
function AbstractIntegrator end


# Print error for integrators not implemented, yet.
function AbstractIntegrator(problem::GeometricProblem, tableau::AbstractTableau; kwargs...)
    error("No integrator found for problem type ", typeof(problem), " and tableau ", tableau)
end


#*****************************************************************************#
# Initialization functions for deterministic integrators                      #
#*****************************************************************************#

# # Create integrator for formal Lagrangian Runge-Kutta tableau.
# function AbstractIntegrator(problem::LODEProblem, tableau::Tableau; kwargs...)
#     IntegratorFLRK(problem, tableau; kwargs...)
# end

# # Create integrator for Projected Gauss-Legendre Runge-Kutta tableau.
# function AbstractIntegrator(problem::Union{IODEProblem,LODEProblem}, tableau::CoefficientsPGLRK; kwargs...)
#     IntegratorPGLRK(problem, tableau; kwargs...)
# end

# # Create integrator for variational partitioned additive Runge-Kutta tableau.
# function AbstractIntegrator(problem::Union{IDAEProblem,LDAEProblem}, tableau::TableauVPARK; kwargs...)
#     IntegratorVPARK(problem, tableau; kwargs...)
# end

# # Create integrator for special partitioned additive Runge-Kutta tableau.
# function AbstractIntegrator(problem::Union{IDAEProblem,LDAEProblem}, tableau::TableauSPARK; kwargs...)
#     IntegratorSPARK(problem, tableau; kwargs...)
# end

# # Create integrator for variational special partitioned additive Runge-Kutta tableau.
# function AbstractIntegrator(problem::Union{IDAEProblem,LDAEProblem}, tableau::TableauVSPARK; kwargs...)
#     IntegratorVSPARK(problem, tableau; kwargs...)
# end

# # Create integrator for variational special partitioned additive Runge-Kutta tableau with projection on primary constraint.
# function AbstractIntegrator(problem::Union{IDAEProblem,LDAEProblem}, tableau::TableauVSPARKprimary; kwargs...)
#     IntegratorVSPARKprimary(problem, tableau; kwargs...)
# end

# # Create integrator for variational special partitioned additive Runge-Kutta tableau with projection on secondary constraint.
# function AbstractIntegrator(problem::LDAEProblem, tableau::TableauVSPARKsecondary; kwargs...)
#     IntegratorVSPARKsecondary(problem, tableau; kwargs...)
# end

# # Create integrator for Hamiltonian partitioned additive Runge-Kutta tableau.
# function AbstractIntegrator(problem::Union{PDAEProblem,HDAEProblem}, tableau::TableauHPARK; kwargs...)
#     IntegratorHPARK(problem, tableau; kwargs...)
# end

# # Create integrator for Hamiltonian special partitioned additive Runge-Kutta tableau.
# function AbstractIntegrator(problem::Union{PDAEProblem,HDAEProblem}, tableau::TableauHSPARK; kwargs...)
#     IntegratorHSPARK(problem, tableau; kwargs...)
# end

# # Create integrator for Hamiltonian special partitioned additive Runge-Kutta tableau with projection on primary constraint.
# function AbstractIntegrator(problem::Union{PDAEProblem,HDAEProblem}, tableau::TableauHSPARKprimary; kwargs...)
#     IntegratorHSPARKprimary(problem, tableau; kwargs...)
# end

# # Create integrator for splitting tableau.
# function AbstractIntegrator(problem::SODEProblem, tableau::AbstractTableauSplitting; kwargs...)
#     IntegratorSplitting(problem, tableau; kwargs...)
# end


#*****************************************************************************#
# Constructor wrappers for deterministic integrators                          #
#*****************************************************************************#

# """
# ```julia
# IntegratorConstructor(DT, D)
# IntegratorConstructor(DT, D, tableau)
# ```

# Returns a closure for various integrator constructors that is needed for
# composition methods based on the tableau type. If not tableau is provided,
# a closure for the constructor of an exact solution is returned.

# - `DT`: data type of the solution / state vector
# - `D`: dimension of the state vector

# """
# function IntegratorConstructor end

# # Create integrator constructor for exact solution.
# function IntegratorConstructor(DT, D)
#     (v::Function, Δt::Number; kwargs...) -> IntegratorExactODE{DT,D}(v, Δt; kwargs...)
# end


#*****************************************************************************#
# General integration functions for all integrators                           #
#*****************************************************************************#

# Apply integrator for ntime time steps and return solution.
function integrate(problem::GeometricProblem, integrator::AbstractIntegrator; kwargs...)
    solution = Solution(problem; kwargs...)
    integrate!(solution, problem, integrator)
    return solution
end

function integrate(problem::GeometricProblem, method::VPRKMethod; kwargs...)
    integrate(problem, AbstractIntegrator(problem, method); kwargs...)
end

# # Integrate given equation with given tableau for ntime time steps and return solution.
# function integrate(problem::GeometricProblem, tableau::AbstractTableau; kwargs...)
#     integrate(problem, AbstractIntegrator(problem, tableau); kwargs...)
# end

# # Integrate ODE specified by vector field and initial condition with given tableau for ntime time steps and return solution.
# function integrate(f::Function, x₀::Vector, tableau::AbstractTableau, tspan, Δt; kwargs...)
#     integrate(ODEProblem(f, tspan, Δt, x₀), tableau; kwargs...)
# end

# # Integrate PODE specified by two vector fields and initial conditions with given tableau for ntime time steps and return solution.
# function integrate(v::Function, f::Function, q₀::Vector, p₀::Vector, tableau::AbstractTableau, tspan, Δt; kwargs...)
#     integrate(PODEProblem(v, f, tspan, Δt, q₀, p₀), tableau; kwargs...)
# end


#*****************************************************************************#
# Integration functions for deterministic integrators                         #
#*****************************************************************************#


# Parts of one integration step that are common to deterministic and stochastic equations.
function integrate!(solstep::SolutionStep, problem::GeometricProblem, int::AbstractIntegrator)
    # integrate one initial condition for one time step
    integrate_step!(int, solstep)

    # take care of periodic solutions
    cut_periodic_solution!(solstep, periodicity(problem))
end

# Integrate equation for initial conditions m with m₁ ≤ m ≤ m₂ for time steps n with n₁ ≤ n ≤ n₂.
function integrate!(sol::GeometricSolution, problem::GeometricProblem, int::AbstractIntegrator, n₁::Int, n₂::Int)
    # check time steps range for consistency
    @assert n₁ ≥ 1
    @assert n₂ ≥ n₁
    @assert n₂ ≤ ntime(sol)

    # create single step solution
    solstep = SolutionStep(problem)

    # copy initial condition from solution
    copy!(solstep, sol[n₁-1])
    initialize!(int, solstep)

    # loop over time steps
    for n in n₁:n₂
        integrate!(solstep, problem, int)
        sol[n] = solstep

        # try
        #     integrate!(int, sol, solstep, n)
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

    return sol
end

# Integrate `equation` for all time steps.
function integrate!(sol::GeometricSolution, problem::GeometricProblem, int::AbstractIntegrator)
    integrate!(sol, problem, int, 1, ntime(sol))
end
