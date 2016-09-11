
"IntegratorIRK: Fully implicit Runge-Kutta integrator."
immutable IntegratorFIRK{T} <: Integrator{T}
    equation::ODE{T}
    tableau::TableauFIRK{T}
    Δt::T

    solver::NonlinearSolver{T}

    x::Array{T,1}
    y::Array{T,1}
    X::Array{T,2}
    Y::Array{T,2}
    F::Array{T,2}

    function IntegratorFIRK(equation, tableau, Δt)
        D = equation.d
        S = tableau.s

        # create solution vectors
        x = zeros(T,D)
        y = zeros(T,D)

        # create internal stage vectors
        X = zeros(T,D,S)
        Y = zeros(T,D,S)
        F = zeros(T,D,S)

        # create solution vector for internal stages / nonlinear solver
        z = zeros(T, D*S)

        # create temporary vectors
        tX = zeros(T,D)
        tF = zeros(T,D)

        # create function parameter datatype


        #
        # tableau.name

        # create function
        code = quote end

        for i in 1:S
            for k in 1:D
                j = D*(i-1)+k
                # copy y to Y
                push!(code.args, :( Y[$k,$i] = y[$j] ))

                # compute X
                push!(code.args, :( X[$k,$i] = x[$k] + params.Δt * Y[$k,$i] ))
            end

            # compute f(X)
            push!(code.args, :( simd_copy_xy_first!(tX, X, $i) ))
            push!(code.args, :( params.f(tX, tF) ))
            push!(code.args, :( simd_copy_yx_first!(tF, F, $i) ))

            # compute b = - (Y-AF)
            for i in 1:S
                for k in 1:D
                    l = D*(i-1)+k
                    push!(code.args, :( b[$l] = - Y[$k,$i] ))
                    for j in 1:size(Y,2)
                        push!(code.args, :( b[$l] += $(tableau.a[i,j]) * F[$k,$j] ))
                    end
                end
            end
        end

        symbolic_name = string(tableau.name)

        # fcode = quote
        @eval begin
            function function_stages!{T, TupleType}(y::Vector{T}, b::Vector{T}, params::NonlinearFunctionParameters{T, Symbol($symbolic_name), TupleType})
                local x::Vector{T} = params.params[1]
                local X::Matrix{T} = params.params[2]
                local Y::Matrix{T} = params.params[3]
                local F::Matrix{T} = params.params[4]
                local tX::Vector{T} = params.params[5]
                local tF::Vector{T} = params.params[6]

                $code
            end
        end

        # println(fcode)
        # eval(fcode)

        params = NonlinearFunctionParameters(tableau.name, equation.f, Δt, x, X, Y, F, tX, tF)

        # create solver
        solver = NewtonSolver(z, params)
        # solver = NewtonSolver(z, equation.f, Δt, x, X, Y, F, tX, tF)
        # TODO allow for other nonlinear solvers based on constructor argument

        # create integrator
        new(equation, tableau, Δt, solver, x, y, X, Y, F)
    end
end

function IntegratorFIRK(equation::Equation, tableau::TableauFIRK, Δt)
    T = eltype(equation.q₀)
    IntegratorFIRK{T}(equation, tableau, Δt)
end

"solve!: Solve ODE with fully implicit Runge-Kutta integrator."
function solve!(int::IntegratorFIRK, sol::SolutionODE)
    local offset::Int

    # copy initial conditions from solution
    simd_copy_xy_first!(int.x, sol, 0)

    for n in 1:sol.ntime
        # compute initial guess
        # TODO
        for i in 1:int.tableau.s
            offset = int.equation.d*(i-1)
            for k in 1:sol.d
                int.solver.x[offset+k] = int.x[k]
            end
        end

        # call nonlinear solver
        solve!(int.solver)
        # println(int.solver.i, ", ", int.solver.rₐ,", ",  int.solver.rᵣ,", ",  int.solver.rₛ)

        # compute final update
        simd_mult!(int.y, int.F, int.tableau.b)
        simd_axpy!(int.Δt, int.y, int.x)

        # copy to solution
        if mod(n, sol.nsave) == 0
            simd_copy_yx_first!(int.x, sol, div(n, sol.nsave))
        end
    end
    nothing
end

"solve!: Solve partitioned ODE with fully implicit Runge-Kutta integrator."
function solve!(int::IntegratorFIRK, s::SolutionPODE)
    # TODO
end
