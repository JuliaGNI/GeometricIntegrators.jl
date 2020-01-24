
#*****************************************************************************#
# General initialization functions for all integrators                        #
#*****************************************************************************#

"Print error for integrators not implemented, yet."
function Integrator(equation::Equation, tableau::AbstractTableau, Δt)
    error("No integrator found for equation ", equation, " and tableau ", tableau)
end


#*****************************************************************************#
# Initialization functions for deterministic integrators                      #
#*****************************************************************************#

"Create integrator for explicit Runge-Kutta tableau."
function Integrator(equation::ODE, tableau::TableauERK, Δt)
    IntegratorERK(equation, tableau, Δt)
end

"Create integrator for diagonally implicit Runge-Kutta tableau."
function Integrator(equation::ODE, tableau::TableauDIRK, Δt)
    IntegratorDIRK(equation, tableau, Δt)
end

"Create integrator for fully implicit Runge-Kutta tableau."
function Integrator(equation::ODE, tableau::TableauFIRK, Δt)
    IntegratorFIRK(equation, tableau, Δt)
end

"Create integrator for explicit partitioned Runge-Kutta tableau."
function Integrator(equation::PODE, tableau::TableauEPRK, Δt)
    IntegratorEPRK(equation, tableau, Δt)
end

"Create integrator for implicit partitioned Runge-Kutta tableau."
function Integrator(equation::PODE, tableau::TableauIPRK, Δt)
    IntegratorIPRK(equation, tableau, Δt)
end

"Create integrator for variational partitioned Runge-Kutta tableau."
function Integrator(equation::IODE, tableau::TableauVPRK, Δt)
    IntegratorVPRK(equation, tableau, Δt)
end

"Create integrator for formal Lagrangian Runge-Kutta tableau."
function Integrator(equation::VODE, tableau::TableauFIRK, Δt)
    IntegratorFLRK(equation, tableau, Δt)
end

"Create integrator for Projected Gauss-Legendre Runge-Kutta tableau."
function Integrator(equation::IODE, tableau::CoefficientsPGLRK, Δt)
    IntegratorPGLRK(equation, tableau, Δt)
end

"Create integrator for variational partitioned additive Runge-Kutta tableau."
function Integrator(equation::IDAE, tableau::TableauVPARK, Δt)
    IntegratorVPARK(equation, tableau, Δt)
end

"Create integrator for variational special partitioned additive Runge-Kutta tableau."
function Integrator(equation::IDAE, tableau::TableauVSPARK, Δt)
    IntegratorVSPARK(equation, tableau, Δt)
end

"Create integrator for variational special partitioned additive Runge-Kutta tableau with projection on primary constraint."
function Integrator(equation::IDAE, tableau::TableauVSPARKprimary, Δt)
    IntegratorVSPARKprimary(equation, tableau, Δt)
end

"Create integrator for variational special partitioned additive Runge-Kutta tableau with projection on secondary constraint."
function Integrator(equation::VDAE, tableau::TableauVSPARKsecondary, Δt)
    IntegratorVSPARKsecondary(equation, tableau, Δt)
end

"Create integrator for Hamiltonian partitioned additive Runge-Kutta tableau."
function Integrator(equation::PDAE, tableau::TableauHPARK, Δt)
    IntegratorHPARK(equation, tableau, Δt)
end

"Create integrator for Hamiltonian special partitioned additive Runge-Kutta tableau."
function Integrator(equation::PDAE, tableau::TableauHSPARK, Δt)
    IntegratorHSPARK(equation, tableau, Δt)
end

"Create integrator for Hamiltonian special partitioned additive Runge-Kutta tableau with projection on primary constraint."
function Integrator(equation::PDAE, tableau::TableauHSPARKprimary, Δt)
    IntegratorHSPARKprimary(equation, tableau, Δt)
end

"Create integrator for splitting tableau."
function Integrator(equation::SODE, tableau::AbstractTableauSplitting, Δt)
    IntegratorSplitting(equation, tableau, Δt)
end


#*****************************************************************************#
# General integration functions for all integrators                           #
#*****************************************************************************#

"Apply integrator for ntime time steps and return solution."
function integrate(integrator::Integrator, ntime::Int, nsave::Int=DEFAULT_NSAVE)
    solution = Solution(equation(integrator), timestep(integrator), ntime, nsave)
    integrate!(integrator, solution)
    return solution
end

"Integrate given equation with given tableau for ntime time steps and return solution."
function integrate(equation::Equation, tableau::AbstractTableau, Δt, ntime, nsave=DEFAULT_NSAVE)
    return integrate(Integrator(equation, tableau, Δt), ntime, nsave)
end

"Integrate ODE specified by vector field and initial condition with given tableau for ntime time steps and return solution."
function integrate(f::Function, x₀::Vector, tableau::AbstractTableau, Δt, ntime, nsave=DEFAULT_NSAVE; t₀=0)
    return integrate(ODE(f, t₀, x₀), tableau, Δt, ntime, nsave)
end

"Integrate PODE specified by two vector fields and initial conditions with given tableau for ntime time steps and return solution."
function integrate(v::Function, f::Function, q₀::Vector, p₀::Vector, tableau::AbstractTableau, Δt, ntime, nsave=DEFAULT_NSAVE; t₀=0)
    return integrate(PODE(v, f, t₀, q₀, p₀), tableau, Δt, ntime, nsave)
end


#*****************************************************************************#
# Integration functions for deterministic integrators                         #
#*****************************************************************************#

"Integrate equation for all initial conditions."
function integrate!(int::DeterministicIntegrator, sol::Solution)
    integrate!(int, sol, 1, sol.ni)
end


"Integrate ODE for initial conditions m with m₁ ≤ m ≤ m₂."
function integrate!(int::DeterministicIntegrator, sol::Solution, m1, m2)
    # integrate initial conditions m with m₁ ≤ m ≤ m₂ for all time steps
    integrate!(int, sol, m1, m2, 1, sol.ntime)
end


"Integrate ODE for initial conditions m with m₁ ≤ m ≤ m₂ for time steps n with n₁ ≤ n ≤ n₂."
function integrate!(int::DeterministicIntegrator{DT,TT}, sol::Solution{DT,TT,N}, m1::Int, m2::Int, n1::Int, n2::Int) where {DT,TT,N}
    @assert m1 ≥ 1
    @assert m2 ≥ m1
    @assert m2 ≤ sol.ni

    @assert n1 ≥ 1
    @assert n2 ≥ n1
    @assert n2 ≤ sol.ntime

    asol = AtomicSolution(equation(int))

    # loop over initial conditions showing progress bar
    for m in m1:m2
        # get cache from solution
        get_initial_conditions!(sol, asol, m, n1)
        initialize!(int, asol)

        # loop over time steps
        for n in n1:n2
            # try
                # integrate one initial condition for one time step
                integrate_step!(int, asol)

                # take care of periodic solutions
                cut_periodic_solution!(asol, periodicity(equation(int)))

                # copy solution from cache to solution
                set_solution!(sol, asol, n, m)
            # catch ex
            #     tstr = " in time step " * string(n)
            #
            #     if m1 ≠ m2
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


# TODO Add solver status information to all integrators (if requested).
