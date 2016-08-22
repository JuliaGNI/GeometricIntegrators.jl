
abstract Integrator

"Integrator: Create integrator for explicit Runge-Kutta tableau."
function Integrator(equation::Equation, tableau::TableauERK)
    IntegratorERK(equation, tableau)
end

"Integrator: Create integrator for diagonally implicit Runge-Kutta tableau."
function Integrator(equation::Equation, tableau::TableauDIRK)
    IntegratorDIRK(equation, tableau)
end

"Integrator: Create integrator for fully implicit Runge-Kutta tableau."
function Integrator(equation::Equation, tableau::TableauFIRK)
    IntegratorFIRK(equation, tableau)
end

"Integrator: Create integrator for partitioned Runge-Kutta tableau."
function Integrator(equation::Equation, tableau::TableauPRK)
    IntegratorPRK(equation, tableau)
end

"Integrator: Create integrator for special additive Runge-Kutta tableau."
function Integrator(equation::Equation, tableau::TableauSARK)
    IntegratorSARK(equation, tableau)
end

"Integrator: Create integrator for special partitioned additive Runge-Kutta tableau."
function Integrator(equation::Equation, tableau::TableauSPARK)
    IntegratorSPARK(equation, tableau)
end

"Integrator: Print error for integrators not implemented, yet."
function Integrator(equation::Equation, tableau::Tableau)
    error("No integrator found for tableau ", tableau)
end

"solve: Solve given equation with given tableau for ntime time steps and return solution."
function solve(equation::Equation, tableau::Tableau, ntime::Int, nsave::Int=1)
    return solve(Integrator(equation, tableau), ntime, nsave)
end

"solve: Apply integrator for ntime time steps and return solution."
function solve(integrator::Integrator, ntime::Int, nsave::Int=1)
    solution = Solution(integrator.equation, ntime, nsave)
    solve!(integrator, solution)
    return solution
end


"IntegratorERK: Explicit Runge-Kutta integrator."
immutable IntegratorERK{T} <: Integrator
    equation::Equation
    tableau::TableauERK

    x::Array{T,1}
    X::Array{T,2}
    Y::Array{T,2}
    F::Array{T,2}

    function IntegratorERK(equation, tableau)
        D = equation.d
        S = tableau.s
        new(equation, tableau, zeros(T,D), zeros(T,D,S), zeros(T,D,S), zeros(T,D,S))
    end
end

function IntegratorERK(equation::Equation, tableau::TableauERK)
    T = eltype(equation.x0)
    IntegratorERK{T}(equation, tableau)
end

"solve!: Solve ODE with explicit Runge-Kutta integrator."
function solve!(int::IntegratorERK, sol::SolutionODE)
    # copy initial conditions from solution
    int.x[:] = sol[1:sol.d, 0]

    for n in 1:sol.ntime
        # compute internal stages
        for i = 1:int.tableau.s
            int.Y[:,i] = 0
            for j = 1:i-1
                int.Y[:,i] += int.tableau.a[i,j] * int.F[:,j]
            end
            int.X[:,i] = int.x[:] + sol.Δt * int.Y[:,i]
            int.F[:,i] = int.equation.f(int.X[:,i])
        end

        # compute final update
        for i in 1:int.tableau.s
            int.x[:] += sol.Δt * int.tableau.b[i] * int.F[:,i]
        end

        # copy to solution
        if mod(n, sol.nsave) == 0
            sol[1:sol.d, div(n, sol.nsave)] = int.x[:]
        end
    end
    return nothing
end

"solve!: Solve partitioned ODE with explicit Runge-Kutta integrator."
function solve!(int::IntegratorERK, s::SolutionPODE)
    # TODO
end


"IntegratorDIRK: Diagonally implicit Runge-Kutta integrator."
immutable IntegratorDIRK{T} <: Integrator
    equation::Equation
    tableau::TableauDIRK

    x::Array{T,1}
    X::Array{T,2}
    Y::Array{T,2}
    F::Array{T,2}

    function IntegratorDIRK(equation, tableau)
        D = equation.d
        S = tableau.s
        new(equation, tableau, zeros(T,D), zeros(T,D,S), zeros(T,D,S), zeros(T,D,S))
    end
end

function IntegratorDIRK(equation::Equation, tableau::TableauDIRK)
    T = eltype(equation.x0)
    IntegratorDIRK{T}(equation, tableau)
end

"solve!: Solve ODE with diagonally implicit Runge-Kutta integrator."
function solve!(int::IntegratorDIRK, s::SolutionODE)
    # TODO
end

"solve!: Solve partitioned ODE with diagonally implicit Runge-Kutta integrator."
function solve!(int::IntegratorDIRK, s::SolutionPODE)
    # TODO
end


"IntegratorIRK: Fully implicit Runge-Kutta integrator."
immutable IntegratorFIRK{T} <: Integrator
    equation::Equation
    tableau::TableauFIRK

    x::Array{T,1}
    X::Array{T,2}
    Y::Array{T,2}
    F::Array{T,2}

    function IntegratorFIRK(equation, tableau)
        D = equation.d
        S = tableau.s
        new(equation, tableau, zeros(T,D), zeros(T,D,S), zeros(T,D,S), zeros(T,D,S))
    end
end

function IntegratorFIRK(equation::Equation, tableau::TableauFIRK)
    T = eltype(equation.x0)
    IntegratorFIRK{T}(equation, tableau)
end

"solve!: Solve ODE with fully implicit Runge-Kutta integrator."
function solve!(int::IntegratorFIRK, s::SolutionODE)
    # TODO
end

"solve!: Solve partitioned ODE with fully implicit Runge-Kutta integrator."
function solve!(int::IntegratorFIRK, s::SolutionPODE)
    # TODO
end


immutable IntegratorPRK{T} <: Integrator
    equation::Equation
    tableau::TableauPRK

end

function solve!(int::IntegratorPRK, s::SolutionPODE)
    # TODO
end


immutable IntegratorSARK{T} <: Integrator
    equation::Equation
    tableau::TableauSARK

end

function solve!(int::IntegratorSARK, s::SolutionDAE)
    # TODO
end


immutable IntegratorSPARK{T} <: Integrator
    equation::Equation
    tableau::TableauSPARK

end

function solve!(int::IntegratorSPARK, s::SolutionPDAE)
    # TODO
end
