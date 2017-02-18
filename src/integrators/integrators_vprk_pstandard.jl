
"Parameters for right-hand side function of variational partitioned Runge-Kutta methods."
type NonlinearFunctionParametersVPRKpStandard{DT,TT,ΑT,GT} <: NonlinearFunctionParameters{DT}
    α::ΑT
    g::GT

    Δt::TT

    d::Int

    R∞::Int

    t::TT

    q̃::Vector{DT}
    p̃::Vector{DT}
    q::Vector{DT}
    p::Vector{DT}
    λ::Vector{DT}

    P::Vector{DT}
    U::Vector{DT}
    G::Vector{DT}

    function NonlinearFunctionParametersVPRKpStandard(α, g, Δt, d, R∞, q̃, p̃)
        @assert d==length(q̃)==length(p̃)

        # create solution vectors
        q = zeros(DT,d)
        p = zeros(DT,d)
        λ = zeros(DT,d)

        P = zeros(DT,d)
        U = zeros(DT,d)
        G = zeros(DT,d)

        new(α, g, Δt, d, R∞, 0, q̃, p̃, q, p, λ, P, U, G)
    end
end

"Compute stages of projected variational partitioned Runge-Kutta methods."
function function_stages!{DT,TT,ΑT,GT}(y::Vector{DT}, b::Vector{DT}, params::NonlinearFunctionParametersVPRKpStandard{DT,TT,ΑT,GT})
    # copy y to q
    for k in 1:params.d
        params.q[k] = y[0*params.d+k]
        params.p[k] = y[1*params.d+k]
        params.λ[k] = y[2*params.d+k]
    end

    # compute p=α(q)
    params.α(params.t+params.Δt, params.q, params.λ, params.P)

    # compute u=λ and g=∇α(q)⋅λ
    simd_copy!(params.λ, params.U)
    params.g(params.t+params.Δt, params.q, params.λ, params.G)

    # compute b = - [q-q̃-U]
    for k in 1:params.d
        b[0*params.d+k] = - params.q[k] + params.q̃[k] + params.R∞ * params.Δt * params.U[k]
    end

    # compute b = - [p-α(q)]
    for k in 1:params.d
        b[1*params.d+k] = - params.p[k] + params.P[k]
    end

    # compute b = - [p-p̃-G]
    for k in 1:params.d
        b[2*params.d+k] = - params.p[k] + params.p̃[k] + params.R∞ * params.Δt * params.G[k]
    end
end

"Variational partitioned Runge-Kutta integrator."
immutable IntegratorVPRKpStandard{DT, TT, ΑT, FT, GT, VT, ST, STP, IT} <: Integrator{DT, TT}
    equation::IODE{DT,TT,ΑT,FT,GT,VT}
    tableau::TableauVPRK{TT}
    Δt::TT

    solver::ST
    projector::STP
    iguess::InitialGuessIODE{DT, TT, VT, FT, IT}

    q::Vector{DT}
    v::Vector{DT}
    p::Vector{DT}
    λ::Vector{DT}

    qₑᵣᵣ::Vector{DT}
    pₑᵣᵣ::Vector{DT}

    y::Vector{DT}
    z::Vector{DT}
    u::Vector{DT}
    g::Vector{DT}

    Q::Array{DT,2}
    V::Array{DT,2}
    U::Array{DT,1}

    P::Array{DT,2}
    F::Array{DT,2}
    G::Array{DT,1}
end

function IntegratorVPRKpStandard{DT,TT,ΑT,FT,GT,VT}(equation::IODE{DT,TT,ΑT,FT,GT,VT}, tableau::TableauVPRK{TT}, Δt::TT;
                                        nonlinear_solver=DEFAULT_NonlinearSolver,
                                        nmax=DEFAULT_nmax, atol=DEFAULT_atol, rtol=DEFAULT_rtol, stol=DEFAULT_stol,
                                        interpolation=HermiteInterpolation{DT})
    D = equation.d
    S = tableau.s

    N = D*S

    if isdefined(tableau, :d)
        d_v = tableau.d
    else
        d_v = DT[]
    end

    # create solver params
    sparams = NonlinearFunctionParametersVPRK{DT,TT,ΑT,FT}(
                                                equation.α, equation.f,
                                                Δt, D, S,
                                                tableau.q, tableau.p,
                                                d_v)

    # create solver
    solver = nonlinear_solver(zeros(DT,N), sparams; nmax=nmax, atol=atol, rtol=rtol, stol=stol)

    # create projector params
    pparams = NonlinearFunctionParametersVPRKpStandard{DT,TT,ΑT,GT}(
                                                equation.α, equation.g,
                                                Δt, D, tableau.R∞, sparams.q, sparams.p)

    # create projector
    projector = nonlinear_solver(zeros(DT,3*D), pparams; nmax=nmax, atol=atol, rtol=rtol, stol=stol)

    # create initial guess
    iguess = InitialGuessIODE(interpolation, equation, Δt)

    # create integrator
    qₑᵣᵣ = zeros(DT,D)
    pₑᵣᵣ = zeros(DT,D)

    y = zeros(DT,D)
    z = zeros(DT,D)
    u = zeros(DT,D)
    g = zeros(DT,D)

    IntegratorVPRKpStandard{DT, TT, ΑT, FT, GT, VT, typeof(solver), typeof(projector), typeof(iguess.int)}(
                                        equation, tableau, Δt, solver, projector, iguess,
                                        sparams.q, sparams.v, sparams.p, pparams.λ,
                                        qₑᵣᵣ, pₑᵣᵣ, y, z, u, g,
                                        sparams.Q, sparams.V, pparams.U,
                                        sparams.P, sparams.F, pparams.G)
end


"Integrate ODE with variational partitioned Runge-Kutta integrator."
function integrate!{DT,TT,ΑT,FT,GT,VT,N}(int::IntegratorVPRKpStandard{DT,TT,ΑT,FT,GT,VT}, sol::SolutionPDAE{DT,TT,N})
    # loop over initial conditions
    for m in 1:sol.ni
        local j::Int
        local tqᵢ::TT
        local tpᵢ::TT

        # copy initial conditions from solution
        get_initial_conditions!(sol, int.q, int.p, m)

        # initialise initial guess
        initialize!(int.iguess, sol.t[0], int.q, int.p)

        for n in 1:sol.ntime
            # set time for nonlinear solver
            int.solver.Fparams.t = sol.t[n]

            # copy previous solution to initial guess
            update!(int.iguess, sol.t[n], int.q, int.p)

            # compute initial guess
            for i in 1:int.tableau.s
                evaluate!(int.iguess, int.y, int.z, int.v, int.tableau.q.c[i], int.tableau.p.c[i])
                for k in 1:int.equation.d
                    int.solver.x[int.equation.d*(i-1)+k] = int.v[k]
                end
            end

            # call nonlinear solver
            solve!(int.solver)

            # println(int.solver.status, ", sit=", n)
            if !solverStatusOK(int.solver.status, int.solver.params)
                println(int.solver.status, ", sit=", n)
            end

            # compute unprojected solution
            simd_mult!(int.y, int.V, int.tableau.q.b)
            simd_mult!(int.z, int.F, int.tableau.p.b)
            simd_axpy!(int.Δt, int.y, int.q, int.qₑᵣᵣ)
            simd_axpy!(int.Δt, int.z, int.p, int.pₑᵣᵣ)

            # set time for projection solver
            int.projector.Fparams.t = sol.t[n]

            # set initial guess for projection
            for k in 1:int.equation.d
                int.projector.x[0*int.equation.d+k] = int.q[k]
                int.projector.x[1*int.equation.d+k] = int.p[k]
                int.projector.x[2*int.equation.d+k] = 0
            end

            # call projection solver
            solve!(int.projector)

            # println(int.projector.status, ", pit=", n)
            if !solverStatusOK(int.projector.status, int.projector.params)
                println(int.projector.status, ", pit=", n)
            end

            if int.solver.status.rₐ == NaN
                break
            end

            # add projection to solution
            simd_copy!(int.U, int.u)
            simd_copy!(int.G, int.g)
            simd_axpy!(int.tableau.R∞*int.Δt, int.u, int.q, int.qₑᵣᵣ)
            simd_axpy!(int.tableau.R∞*int.Δt, int.g, int.p, int.pₑᵣᵣ)

            # copy to solution
            copy_solution!(sol, int.q, int.p, int.λ, n, m)
        end
    end
    nothing
end
