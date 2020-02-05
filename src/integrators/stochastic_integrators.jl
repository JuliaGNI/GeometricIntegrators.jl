
#*****************************************************************************#
# Initialization functions for stochastic integrators                         #
#*****************************************************************************#

"Create integrator for stochastic explicit Runge-Kutta tableau."
function Integrator(equation::SDE, tableau::TableauSERK, Δt)
    IntegratorSERK(equation, tableau, Δt)
end

"Create integrator for weak explicit Runge-Kutta tableau."
function Integrator(equation::SDE, tableau::TableauWERK, Δt)
    IntegratorWERK(equation, tableau, Δt)
end

"Create integrator for stochastic fully implicit Runge-Kutta tableau."
function Integrator(equation::SDE, tableau::TableauSIRK, Δt; K::Int=0)
    IntegratorSIRK(equation, tableau, Δt, K=K)
end

"Create integrator for stochastic fully implicit partitioned Runge-Kutta tableau."
function Integrator(equation::PSDE, tableau::TableauSIPRK, Δt; K::Int=0)
    IntegratorSIPRK(equation, tableau, Δt, K=K)
end

"Create integrator for stochastic fully implicit split partitioned Runge-Kutta tableau."
function Integrator(equation::SPSDE, tableau::TableauSISPRK, Δt; K::Int=0)
    IntegratorSISPRK(equation, tableau, Δt, K=K)
end

"Create integrator for weak fully implicit Runge-Kutta tableau."
function Integrator(equation::SDE, tableau::TableauWIRK, Δt)
    IntegratorWIRK(equation, tableau, Δt)
end


#*****************************************************************************#
# Integration functions for stochastic integrators                            #
#*****************************************************************************#

"Integrate SDE for all sample paths and initial conditions."
function integrate!(int::StochasticIntegrator, sol::StochasticSolution)
    integrate!(int, sol, 1, nsamples(sol))
end


"Integrate SDE for the sample paths m with m₁ ≤ m ≤ m₂."
function integrate!(int::StochasticIntegrator, sol::StochasticSolution, m1::Int, m2::Int)
    # integrate sample paths m with m₁ ≤ m ≤ m₂ for all time steps
    integrate!(int, sol, m1, m2, 1, ntime(sol))
end


"Integrate SDE for the sample paths m with m₁ ≤ m ≤ m₂, and the time steps n with n₁ ≤ n ≤ n₂."
function integrate!(int::StochasticIntegrator{DT,TT}, sol::StochasticSolution{DT,TT,N}, m1::Int, m2::Int, n1::Int, n2::Int) where {DT,TT,N}
    @assert m1 ≥ 1
    @assert m2 ≥ m1
    @assert m2 ≤ nsamples(sol)

    @assert n1 ≥ 1
    @assert n2 ≥ n1
    @assert n2 ≤ ntime(sol)

    asol = AtomicSolution(equation(int))

    # loop over sample paths
    for m in m1:m2
        get_initial_conditions!(sol, asol, m, n1)
        initialize!(int, asol)

        # loop over time steps
        for n in n1:n2
            # try
            integrate!(int, asol, m, n)
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
            #         throw(ex)
            #     else
            #         @warn(string(typeof(ex)) * tstr)
            #         throw(ex)
            #     end
            # end
        end
    end
end


function integrate!(int::StochasticIntegrator{DT,TT}, asol::AtomicSolution{DT,TT}, m::Int, n::Int) where {DT,TT}
    # copy the increments of the Brownian Process
    get_increments!(sol, asol, n, m)

    # integrate one initial condition for one time step
    integrate_step!(int, asol)

    # take care of periodic solutions
    cut_periodic_solution!(asol, periodicity(equation(int)))

    # copy solution from cache to solution
    set_solution!(sol, asol, n, m)
end
