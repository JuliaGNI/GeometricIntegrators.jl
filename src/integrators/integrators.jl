
abstract Integrator

"Integrator: Create integrator for explicit Runge-Kutta tableau."
function Integrator(equation::Equation, tableau::TableauERK, Δt)
    IntegratorERK(equation, tableau, Δt)
end

"Integrator: Create integrator for diagonally implicit Runge-Kutta tableau."
function Integrator(equation::Equation, tableau::TableauDIRK, Δt)
    IntegratorDIRK(equation, tableau, Δt)
end

"Integrator: Create integrator for fully implicit Runge-Kutta tableau."
function Integrator(equation::Equation, tableau::TableauFIRK, Δt)
    IntegratorFIRK(equation, tableau, Δt)
end

"Integrator: Create integrator for partitioned Runge-Kutta tableau."
function Integrator(equation::Equation, tableau::TableauPRK, Δt)
    IntegratorPRK(equation, tableau, Δt)
end

"Integrator: Create integrator for special additive Runge-Kutta tableau."
function Integrator(equation::Equation, tableau::TableauSARK, Δt)
    IntegratorSARK(equation, tableau, Δt)
end

"Integrator: Create integrator for special partitioned additive Runge-Kutta tableau."
function Integrator(equation::Equation, tableau::TableauSPARK, Δt)
    IntegratorSPARK(equation, tableau, Δt)
end

"Integrator: Print error for integrators not implemented, yet."
function Integrator(equation::Equation, tableau::Tableau, Δt)
    error("No integrator found for tableau ", tableau)
end

"solve: Solve given equation with given tableau for ntime time steps and return solution."
function solve(equation::Equation, tableau::Tableau, Δt, ntime::Int, nsave::Int=1)
    return solve(Integrator(equation, tableau, Δt), ntime, nsave)
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
    Δt::T

    x::Array{T,1}
    X::Array{T,2}
    Y::Array{T,2}
    F::Array{T,2}

    function IntegratorERK(equation, tableau, Δt)
        D = equation.d
        S = tableau.s
        new(equation, tableau, Δt, zeros(T,D), zeros(T,D,S), zeros(T,D,S), zeros(T,D,S))
    end
end

function IntegratorERK(equation::Equation, tableau::TableauERK, Δt)
    T = eltype(equation.q₀)
    IntegratorERK{T}(equation, tableau, Δt)
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
            int.X[:,i] = int.x[:] + int.Δt * int.Y[:,i]
            # int.equation.f(int.X[:,i], int.F[:,i])
            int.equation.f(view(int.X, 1:sol.d, i), view(int.F, 1:sol.d, i))
        end

        # compute final update
        for i in 1:int.tableau.s
            int.x[:] += int.Δt * int.tableau.b[i] * int.F[:,i]
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
    Δt::T

    x::Array{T,1}
    X::Array{T,2}
    Y::Array{T,2}
    F::Array{T,2}

    function IntegratorDIRK(equation, tableau, Δt)
        D = equation.d
        S = tableau.s
        new(equation, tableau, Δt, zeros(T,D), zeros(T,D,S), zeros(T,D,S), zeros(T,D,S))
    end
end

function IntegratorDIRK(equation::Equation, tableau::TableauDIRK, Δt)
    T = eltype(equation.q₀)
    IntegratorDIRK{T}(equation, tableau, Δt)
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
    Δt::T

    solver::NonlinearSolver{T}

    x::Array{T,1}
    X::Array{T,2}
    Y::Array{T,2}
    F::Array{T,2}

    function IntegratorFIRK(equation, tableau, Δt)
        D = equation.d
        S = tableau.s

        # create solution vector
        x = zeros(T,D)

        # create internal stage vectors
        X = zeros(T,D,S)
        Y = zeros(T,D,S)
        F = zeros(T,D,S)

        # create solution vector for internal stages / nonlinear solver
        z = zeros(T, D*S)

        # create function
        function function_stages{T}(y::Vector{T}, b::Vector{T})
            # loop through stages
            for i in 1:S
                # copy y to Y
                Y[:,i] = y[D*(i-1)+1:D*i]

                # compute X
                X[:,i] = x[:] + Δt * Y[:,i]

                # compute f(X)
                # equation.f(X[:,i], F[:,i])
                equation.f(view(X, 1:D, i), view(F, 1:D, i))
            end

            # compute b = - (Y-AF)
            for i in 1:S
                b[D*(i-1)+1:D*i] = - Y[:,i]
                for j in 1:S
                    b[D*(i-1)+1:D*i] += tableau.a[i,j] * F[:,j]
                end
            end
        end

        # create solver
        solver = NewtonSolver(z, function_stages)

        # create integrator
        new(equation, tableau, Δt, solver, x, X, Y, F)
    end
end

function IntegratorFIRK(equation::Equation, tableau::TableauFIRK, Δt)
    T = eltype(equation.q₀)
    IntegratorFIRK{T}(equation, tableau, Δt)
end

"solve!: Solve ODE with fully implicit Runge-Kutta integrator."
function solve!(int::IntegratorFIRK, sol::SolutionODE)
    # copy initial conditions from solution
    int.x[:] = sol[1:sol.d, 0]

    for n in 1:sol.ntime
        # compute initial guess
        # TODO
        for i in 1:int.tableau.s
            int.solver.z[int.equation.d*(i-1)+1:int.equation.d*i] = int.x[:]
        end

        # call nonlinear solver
        solve!(int.solver)
        # println(int.solver.i, ", ", int.solver.rₐ,", ",  int.solver.rᵣ,", ",  int.solver.rₛ)

        # compute final update
        for i in 1:int.tableau.s
            int.x[:] += int.Δt * int.tableau.b[i] * int.F[:,i]
        end

        # copy to solution
        if mod(n, sol.nsave) == 0
            sol[1:sol.d, div(n, sol.nsave)] = int.x[:]
        end
    end
    return nothing
end

# "solve!: Solve partitioned ODE with fully implicit Runge-Kutta integrator."
# function solve!(int::IntegratorFIRK, s::SolutionPODE)
#     # TODO
# end


"IntegratorPRK: Explicit partitioned Runge-Kutta integrator."
immutable IntegratorPRK{T} <: Integrator
    equation::Equation
    tableau::TableauPRK
    Δt::T

    q::Array{T,1}
    p::Array{T,1}
    Q::Array{T,2}
    P::Array{T,2}
    Y::Array{T,2}
    Z::Array{T,2}
    F::Array{T,2}
    G::Array{T,2}

    function IntegratorPRK(equation, tableau, Δt)
        D = equation.d
        S = tableau.s
        new(equation, tableau, Δt,
            zeros(T,D), zeros(T,D),
            zeros(T,D,S), zeros(T,D,S),
            zeros(T,D,S), zeros(T,D,S),
            zeros(T,D,S), zeros(T,D,S))
    end
end

function IntegratorPRK(equation::Equation, tableau::TableauPRK, Δt)
    T1 = eltype(equation.q₀)
    T2 = eltype(equation.p₀)
    @assert T1 == T2
    IntegratorPRK{T1}(equation, tableau, Δt)
end


function computeStageQ!(int::IntegratorPRK, Δt::AbstractFloat, i::Integer, jmax::Integer)
    d = length(int.Y[:,i])
    int.Y[:,i] = 0
    for j = 1:jmax
        int.Y[:,i] += int.tableau.a_q[i,j] * int.F[:,j]
    end
    int.Q[:,i] = int.q[:] + Δt * int.Y[:,i]
    # int.equation.g(int.Q[:,i], int.G[:,i])
    int.equation.g(view(int.Q, 1:d, i), view(int.G, 1:d, i))
end

function computeStageP!(int::IntegratorPRK, Δt::AbstractFloat, i::Integer, jmax::Integer)
    d = length(int.Z[:,i])
    int.Z[:,i] = 0
    for j = 1:jmax
        int.Z[:,i] += int.tableau.a_p[i,j] * int.G[:,j]
    end
    int.P[:,i] = int.p[:] + Δt * int.Z[:,i]
    # int.equation.f(int.P[:,i], int.F[:,i])
    int.equation.f(view(int.P, 1:d, i), view(int.F, 1:d, i))
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
                computeStageP!(int, int.Δt, i, i-1)
                computeStageQ!(int, int.Δt, i, i)
            elseif int.tableau.a_p[i,i] ≠ 0.
                computeStageQ!(int, int.Δt, i, i-1)
                computeStageP!(int, int.Δt, i, i)
            else
                computeStageQ!(int, int.Δt, i, i-1)
                computeStageP!(int, int.Δt, i, i-1)
            end
        end

        # compute final update
        for i in 1:int.tableau.s
            int.q[:] += int.Δt * int.tableau.b_q[i] * int.F[:,i]
            int.p[:] += int.Δt * int.tableau.b_p[i] * int.G[:,i]
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
