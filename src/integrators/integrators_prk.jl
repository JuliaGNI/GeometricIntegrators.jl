
"IntegratorPRK: Explicit partitioned Runge-Kutta integrator."
immutable IntegratorPRK{T} <: Integrator{T}
    equation::PODE{T}
    tableau::TableauPRK{T}
    Δt::T

    q::Array{T,1}
    p::Array{T,1}
    y::Array{T,1}
    z::Array{T,1}
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


function computeStageQ!(int::IntegratorPRK, i::Int, jmax::Int)
    for j in 1:jmax
        for k in 1:int.equation.d
            int.Y[k,i] += int.tableau.a_q[i,j] * int.F[k,j]
        end
    end
    for k in 1:int.equation.d
        int.Q[k,i] = int.q[k] + int.Δt * int.Y[k,i]
    end
    int.equation.g(view(int.Q, :, i), view(int.G, :, i))
    nothing
end

function computeStageP!(int::IntegratorPRK, i::Int, jmax::Int)
    for j in 1:jmax
        for k in 1:int.equation.d
            int.Z[k,i] += int.tableau.a_p[i,j] * int.G[k,j]
        end
    end
    for k in 1:int.equation.d
        int.P[k,i] = int.p[k] + int.Δt * int.Z[k,i]
    end
    int.equation.f(view(int.P, :, i), view(int.F, :, i))
    nothing
end

"solve!: Solve partitioned ODE with explicit partitioned Runge-Kutta integrator."
function solve!(int::IntegratorPRK, sol::SolutionPODE)
    local j::Int
    # copy initial conditions from solution
    int.q .= sol[1:sol.d, 1, 0]
    int.p .= sol[1:sol.d, 2, 0]
    stages = 1:int.tableau.s

    for n in 1:sol.ntime
        # compute internal stages
        fill!(int.Y, 0.)
        fill!(int.Z, 0.)
        for i in stages
            if int.tableau.a_q[i,i] ≠ 0. && int.tableau.a_p[i,i] ≠ 0.
                error("This is an implicit method!")
            elseif int.tableau.a_q[i,i] ≠ 0.
                computeStageP!(int, i, i-1)
                computeStageQ!(int, i, i)
            elseif int.tableau.a_p[i,i] ≠ 0.
                computeStageQ!(int, i, i-1)
                computeStageP!(int, i, i)
            else
                computeStageQ!(int, i, i-1)
                computeStageP!(int, i, i-1)
            end
        end

        # compute final update
        simd_mult!(int.y, int.F, int.tableau.b_q)
        simd_axpy!(int.Δt, int.y, int.q)

        simd_mult!(int.z, int.G, int.tableau.b_p)
        simd_axpy!(int.Δt, int.z, int.p)

        # copy to solution
        if mod(n, sol.nsave) == 0
            j = div(n, sol.nsave)
            for i in 1:sol.d
                sol[i, 1, j] = int.q[i]
                sol[i, 2, j] = int.p[i]
            end
        end
    end
    nothing
end
