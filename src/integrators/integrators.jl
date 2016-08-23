
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


# TODO Add solver status information to all integrators.


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


"IntegratorPRK: Explicit partitioned Runge-Kutta integrator."
immutable IntegratorPRK{T} <: Integrator
    equation::Equation
    tableau::TableauPRK

    q::Array{T,1}
    p::Array{T,1}
    Q::Array{T,2}
    P::Array{T,2}
    Y::Array{T,2}
    Z::Array{T,2}
    F::Array{T,2}
    G::Array{T,2}

    function IntegratorPRK(equation, tableau)
        D = equation.d
        S = tableau.s
        new(equation, tableau,
            zeros(T,D), zeros(T,D),
            zeros(T,D,S), zeros(T,D,S),
            zeros(T,D,S), zeros(T,D,S),
            zeros(T,D,S), zeros(T,D,S))
    end
end

function IntegratorPRK(equation::Equation, tableau::TableauPRK)
    T1 = eltype(equation.q0)
    T2 = eltype(equation.p0)
    @assert T1 == T2
    IntegratorPRK{T1}(equation, tableau)
end


function computeStageQ!(int::IntegratorPRK, Δt::AbstractFloat, i::Integer, jmax::Integer)
    int.Y[:,i] = 0
    for j = 1:jmax
        int.Y[:,i] += int.tableau.a_q[i,j] * int.F[:,j]
    end
    int.Q[:,i] = int.q[:] + Δt * int.Y[:,i]
    int.G[:,i] = int.equation.g(int.Q[:,i])
end

function computeStageP!(int::IntegratorPRK, Δt::AbstractFloat, i::Integer, jmax::Integer)
    int.Z[:,i] = 0
    for j = 1:jmax
        int.Z[:,i] += int.tableau.a_p[i,j] * int.G[:,j]
    end
    int.P[:,i] = int.p[:] + Δt * int.Z[:,i]
    int.F[:,i] = int.equation.f(int.P[:,i])
end

"solve!: Solve partitioned ODE with explicit partitioned Runge-Kutta integrator."
function solve!(int::IntegratorPRK, sol::SolutionPODE)
    # copy initial conditions from solution
    int.q[:] = sol[1:sol.d, 1, 0]
    int.p[:] = sol[1:sol.d, 2, 0]

    for n in 1:sol.ntime
        # compute internal stages
        for i = 1:int.tableau.s
            if int.tableau.a_q[i,i] ≠ 0. && int.tableau.a_p[i,i] ≠ 0.
                error("This is an implicit method!")
            elseif int.tableau.a_q[i,i] ≠ 0.
                computeStageP!(int, sol.Δt, i, i-1)
                computeStageQ!(int, sol.Δt, i, i)
            elseif int.tableau.a_p[i,i] ≠ 0.
                computeStageQ!(int, sol.Δt, i, i-1)
                computeStageP!(int, sol.Δt, i, i)
            else
                computeStageQ!(int, sol.Δt, i, i-1)
                computeStageP!(int, sol.Δt, i, i-1)
            end
        end

        # compute final update
        for i in 1:int.tableau.s
            int.q[:] += sol.Δt * int.tableau.b_q[i] * int.F[:,i]
            int.p[:] += sol.Δt * int.tableau.b_p[i] * int.G[:,i]
        end

        # copy to solution
        if mod(n, sol.nsave) == 0
            sol[1:sol.d, 1, div(n, sol.nsave)] = int.q[:]
            sol[1:sol.d, 2, div(n, sol.nsave)] = int.p[:]
        end
    end
    return nothing
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
