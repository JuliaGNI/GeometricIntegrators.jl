
abstract type Integrator{dType, tType} end


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

"Create integrator for singly implicit Runge-Kutta tableau."
function Integrator(equation::ODE, tableau::TableauSIRK, Δt)
    IntegratorSIRK(equation, tableau, Δt)
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

"Create integrator for additive Runge-Kutta tableau."
function Integrator(equation::DAE, tableau::TableauARK, Δt)
    IntegratorARK(equation, tableau, Δt)
end

"Create integrator for special additive Runge-Kutta tableau."
function Integrator(equation::DAE, tableau::TableauSARK, Δt)
    IntegratorSARK(equation, tableau, Δt)
end

"Create integrator for partitioned additive Runge-Kutta tableau."
function Integrator(equation::PDAE, tableau::TableauPARK, Δt)
    IntegratorPARK(equation, tableau, Δt)
end

"Create integrator for special partitioned additive Runge-Kutta tableau."
function Integrator(equation::PDAE, tableau::TableauSPARK, Δt)
    IntegratorSPARK(equation, tableau, Δt)
end

"Create integrator for variational partitioned additive Runge-Kutta tableau."
function Integrator(equation::IDAE, tableau::TableauVPARK, Δt)
    IntegratorVPARK(equation, tableau, Δt)
end

"Create integrator for variational special partitioned additive Runge-Kutta tableau."
function Integrator(equation::IDAE, tableau::TableauVSPARK, Δt)
    IntegratorVSPARK(equation, tableau, Δt)
end

"Print error for integrators not implemented, yet."
function Integrator(equation::Equation, tableau::AbstractTableau, Δt)
    error("No integrator found for tableau ", tableau)
end

"Apply integrator for ntime time steps and return solution."
function integrate(integrator::Integrator, ntime::Int, nsave::Int=1)
    solution = Solution(integrator.equation, integrator.Δt, ntime, nsave)
    integrate!(integrator, solution)
    return solution
end

"Integrate given equation with given tableau for ntime time steps and return solution."
function integrate(equation::Equation, tableau::AbstractTableau, Δt, ntime::Int, nsave::Int=1)
    return integrate(Integrator(equation, tableau, Δt), ntime, nsave)
end

"Integrate ODE specified by vector field and initial condition with given tableau for ntime time steps and return solution."
function integrate(f::Function, x₀::Vector, tableau::AbstractTableau, Δt, ntime::Int, nsave::Int=1; t₀=0)
    return integrate(ODE(f, t₀, x₀), tableau, Δt, ntime, nsave)
end


"Initialize integrator for initial conditions m with m₁ ≤ m ≤ m₂ and time step 0."
function initialize!(int::Integrator, sol::Solution, m1::Int, m2::Int)
    for m in m1:m2
        # initialize integrator for initial condition m and time step 0
        initialize!(int, sol, m)
    end
end

"Integrate ODE for all initial conditions."
function integrate!(int, sol)
    integrate!(int, sol, 1, sol.ni)
end

# TODO Add counter to solution and reactivate this.
# "Integrate ODE for all initial conditions for nt time steps."
# function integrate!(int, sol, ntime)
#     integrate!(int, sol, 1, sol.ni, ntime)
# end

"Integrate ODE for initial conditions m with m₁ ≤ m ≤ m₂."
function integrate!(int, sol, m1, m2)
    # initialize integrator for initial conditions m with m₁ ≤ m ≤ m₂ and time step 0
    initialize!(int, sol, m1, m2)

    # integrate initial conditions m with m₁ ≤ m ≤ m₂ for all time steps
    integrate!(int, sol, m1, m2, 1, sol.ntime)
end

"Integrate ODE for initial conditions m with m₁ ≤ m ≤ m₂ for time steps n with n₁ ≤ n ≤ n₂."
function integrate!(int::Integrator{DT,TT}, sol::Solution{DT,TT,N}, m1::Int, m2::Int, n1::Int, n2::Int) where {DT,TT,N}
    @assert m1 ≥ 1
    @assert m2 ≥ m1
    @assert m2 ≤ sol.ni

    @assert n1 ≥ 1
    @assert n2 ≥ n1
    @assert n2 ≤ sol.ntime

    local nshow = get_config(:int_show_progress_nmin)
    local nrun  = (m2-m1+1)*(n2-n1+1)

    # initialize progress bar
    if nrun ≥ nshow
        p = Progress(nrun, 5)
    end

    # loop over initial conditions
    for m in m1:m2
        # loop over time steps
        for n in n1:n2
            # integrate one initial condition for one time step
            integrate_step!(int, sol, m, n)

            # update progress bar
            if nrun ≥ nshow
                next!(p)
            end
        end
    end
end


# TODO Add solver status information to all integrators (if requested).
