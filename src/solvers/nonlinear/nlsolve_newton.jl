
using NLsolve: OnceDifferentiable, NewtonCache, newton_
using LineSearches


struct NLsolveNewton{T, FT, DT, CT, ST, LT} <: AbstractNewtonSolver{T}
    x::Vector{T}
    f::Vector{T}
    J::Matrix{T}

    cache::CT

    F!::FT
    DF::DT

    line_search::ST
    linear_solver::LT

    params::NonlinearSolverParameters{T}
    status::NonlinearSolverStatus{T}

    function NLsolveNewton(x::AbstractVector{T}, f::AbstractVector{T}, J::AbstractMatrix{T},
                    F!::FT, DF::DT, cache::CT, line_search::ST, linear_solver::LT) where {T,FT,DT,CT,ST,LT}

        nls_params = NonlinearSolverParameters(T)
        nls_status = NonlinearSolverStatus{T}(length(x))

        new{T,FT,DT,CT,ST,LT}(x, f, J, cache, F!, DF, line_search, linear_solver, nls_params, nls_status)
    end
end


function NLsolveNewton(x::AbstractVector{T}, F!::Function; J!::Union{Function,Nothing}=nothing) where {T}
    linear_solver = getLinearSolver(x)

    f  = similar(x)
    df = zero(linear_solver.A)

    f! = (f, x) -> F!(x, f)

    if J! == nothing
        if get_config(:jacobian_autodiff)
            DF = OnceDifferentiable(f!, x, f, df, :forward)
        else
            DF = OnceDifferentiable(f!, x, f, df, :finite)
        end
    else
        j! = (df, x) -> J!(x, df)
        DF = OnceDifferentiable(f!, j!, x, f, df)
    end

    NLsolveNewton(x, f, df, F!, DF, NewtonCache(DF), LineSearches.Static(), linear_solver)
end


function solve!(s::NLsolveNewton; n::Int=0)
    linsolve = (x, A, b) -> begin
        s.linear_solver.A .= A
        s.linear_solver.b .= b
        factorize!(s.linear_solver)
        solve!(s.linear_solver)
        x .= s.linear_solver.b
    end

    res=newton_(s.DF, s.x, s.params.stol, s.params.atol, s.params.nmax, false, false, false,
                s.line_search, linsolve, s.cache)

    s.x .= res.zero
    s.status.i  = res.iterations
    s.status.rₐ = res.residual_norm

    nothing
end
