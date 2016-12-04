
using GeomDAE.Fields: Field, Q, P, V, F

"Parameters for right-hand side function of implicit partitioned Runge-Kutta methods."
type NonlinearFunctionParametersIPRK{T} <: NonlinearFunctionParameters{T}
    fields::Tuple{Vararg{Field}}
    functs::Tuple{Vararg{Field}}

    fq::Function
    fp::Function
    fv::Function
    ff::Function

    Δt::T
    t::T

    d::Int
    s::Int
    a_q::Matrix{T}
    a_p::Matrix{T}
    c_q::Vector{T}
    c_p::Vector{T}

    q::Vector{T}
    p::Vector{T}
    y::Vector{T}
    z::Vector{T}

    Q::Matrix{T}
    V::Matrix{T}
    P::Matrix{T}
    F::Matrix{T}

    Y::Array{T,3}

    tQ::Vector{T}
    tV::Vector{T}
    tP::Vector{T}
    tF::Vector{T}

    function NonlinearFunctionParametersIPRK(fields, functs, fq, fp, fv, ff, Δt, d, s, a_q, a_p, c_q, c_p)
        # create solution vectors
        q = zeros(T,d)
        p = zeros(T,d)
        y = zeros(T,d)
        z = zeros(T,d)

        # create internal stage vectors
        Q = zeros(T,d,s)
        V = zeros(T,d,s)
        P = zeros(T,d,s)
        F = zeros(T,d,s)

        # create temporary vectors
        Y  = zeros(T,2,d,s)

        tQ = zeros(T,d)
        tV = zeros(T,d)
        tP = zeros(T,d)
        tF = zeros(T,d)

        new(fields, functs, fq, fp, fv, ff, Δt, 0, d, s, a_q, a_p, c_q, c_p, q, p, y, z, Q, V, P, F, Y, tQ, tV, tP, tF)
    end
end

"Compute stages of implicit partitioned Runge-Kutta methods."
@generated function function_stages!{T}(y::Vector{T}, b::Vector{T},
                             params::NonlinearFunctionParametersIPRK{T},
                             fields::Field...)

    local fargs = []

     if Q ∈ fields
         push!(fargs, :(params.tQ))
     end

     if V ∈ fields
         push!(fargs, :(params.tV))
     end

     if P ∈ fields
         push!(fargs, :(params.tP))
     end

     if F ∈ fields
         push!(fargs, :(params.tF))
     end

     # TODO Below, only the case (fp, ff) is implemented. All other combinations
     #      need to be implemented as well.

     code = quote
         local tqᵢ::T
         local tpᵢ::T

         # loop through stages
         for i in 1:params.s
             for k in 1:params.d
                 # copy y to Y
                 params.Y[1,k,i] = y[2*(params.d*(i-1)+k-1)+1]
                 params.Y[2,k,i] = y[2*(params.d*(i-1)+k-1)+2]

                 # compute Q and V
                 params.Q[k,i] = params.q[k] + params.Δt * params.Y[1,k,i]
                 params.V[k,i] = params.Y[2,k,i]
             end

             # compute f(X)
             tqᵢ = params.t + params.Δt * params.c_q[i]
             tpᵢ = params.t + params.Δt * params.c_p[i]
             simd_copy_xy_first!(params.tQ, params.Q, i)
             simd_copy_xy_first!(params.tV, params.V, i)
             params.fp(tpᵢ, params.tQ, params.tV, params.tP)
             params.ff(tpᵢ, params.tQ, params.tV, params.tF)
             simd_copy_yx_first!(params.tP, params.P, i)
             simd_copy_yx_first!(params.tF, params.F, i)
         end

         # compute b = - [(Y-AV), (P-AF)]
         for i in 1:params.s
             for k in 1:params.d
                 b[2*(params.d*(i-1)+k-1)+1] = - params.Y[1,k,i]
                 b[2*(params.d*(i-1)+k-1)+2] = - params.P[k,i] + params.p[k]
                 for j in 1:params.s
                     b[2*(params.d*(i-1)+k-1)+1] += params.a_q[i,j] * params.V[k,j]
                     b[2*(params.d*(i-1)+k-1)+2] += params.a_p[i,j] * params.F[k,j] * params.Δt
                 end
             end
         end
    end

    return code
end

"Compute stages of implicit partitioned Runge-Kutta methods."
function function_stages!{T}(y::Vector{T}, b::Vector{T},
                             params::NonlinearFunctionParametersIPRK{T})
    function_stages!(y, b, params, params.fields...)
end


"Implicit partitioned Runge-Kutta integrator."
immutable IntegratorIPRK{T, ST} <: Integrator{T}
    equation::SODE{T}
    tableau::TableauIPRK{T}
    Δt::T

    solver::ST

    q::Array{T,1}
    p::Array{T,1}
    y::Array{T,1}
    z::Array{T,1}
    Q::Array{T,2}
    V::Array{T,2}
    P::Array{T,2}
    F::Array{T,2}
end

function IntegratorIPRK{T}(equation::SODE{T}, tableau::TableauIPRK{T}, Δt::T; nonlinear_solver=QuasiNewtonSolver)
    D = equation.d
    S = tableau.s

    # create solution vector for internal stages / nonlinear solver
    z = zeros(T, 2*D*S)

    # create params
    params = NonlinearFunctionParametersIPRK{T}(equation.fields, equation.functs,
                                                equation.q, equation.p, equation.v, equation.f,
                                                Δt, D, S,
                                                tableau.a_q, tableau.a_p, tableau.c_q, tableau.c_p)

    # create solver
    solver = nonlinear_solver(z, params)

    # create integrator
    IntegratorIPRK{T, typeof(solver)}(equation, tableau, Δt, solver,
                                      params.q, params.p, params.y, params.z,
                                      params.Q, params.V, params.P, params.F)
end


"Integrate ODE with special partitioned Runge-Kutta integrator."
function integrate!(int::IntegratorIPRK, sol::SolutionPODE)
    local nt::Int

    # copy initial conditions from solution
    simd_copy_xy_first!(int.q, sol.q, 1)
    simd_copy_xy_first!(int.p, sol.p, 1)

    for n in 1:sol.ntime
        # compute initial guess
        # TODO
        for i in 1:int.tableau.s
            for k in 1:int.equation.d
                int.solver.x[2*(sol.nd*(i-1)+k-1)+1] = int.q[k]
                # TODO initial guess for velocity
                # int.solver.x[2*(params.d*(i-1)+k-1)+1] = 0.
            end
        end

        # call nonlinear solver
        solve!(int.solver)

        if !solverStatusOK(int.solver.status, int.solver.params)
            println(int.solver.status)
        end

        # compute final update
        simd_mult!(int.y, int.V, int.tableau.b_q)
        simd_mult!(int.z, int.F, int.tableau.b_p)
        simd_axpy!(int.Δt, int.y, int.q)
        simd_axpy!(int.Δt, int.z, int.p)

        # copy to solution
        if mod(n, sol.nsave) == 0
            nt = div(n, sol.nsave)+1
            simd_copy_yx_first!(int.q, sol.q, nt)
            simd_copy_yx_first!(int.p, sol.p, nt)
        end
    end
    nothing
end
