
"Parameters for right-hand side function of Specialised Partitioned Additive Runge-Kutta methods for Variational systems."
mutable struct AbstractParametersVSPARK{IT,DT,TT,D,S,R,P,FT,PT,UT,GT,ϕT,tabType} <: AbstractParametersSPARK{DT,TT}
    f_f::FT
    f_p::PT
    f_u::UT
    f_g::GT
    f_ϕ::ϕT

    Δt::TT

    tab::tabType

    @ParametersSPARK

    function AbstractParametersVSPARK{IT,DT,D,S,R,P}(f_f::FT, f_p::PT, f_u::UT, f_g::GT, f_ϕ::ϕT, Δt::TT, tableau::tabType) where {IT,DT,TT,D,S,R,P,FT,PT,UT,GT,ϕT,tabType}
        # create solution vectors
        q = zeros(DT,D)
        p = zeros(DT,D)
        λ = zeros(DT,D)

        new{IT,DT,TT,D,S,R,P,FT,PT,UT,GT,ϕT,tabType}(f_f, f_p, f_u, f_g, f_ϕ, Δt, tableau, zero(TT), q, p, λ)
    end
end


@inline equation(int::AbstractIntegratorVSPARK) = int.equation
@inline timestep(int::AbstractIntegratorVSPARK) = int.params.Δt
@inline tableau(int::AbstractIntegratorVSPARK) = int.tableau
@inline pstages(int::AbstractIntegratorVSPARK) = int.tableau.r


function initial_guess!(int::AbstractIntegratorVSPARK, sol::AtomicSolutionPDAE)
    for i in eachstage(int)
        evaluate!(int.iguess, sol.q, sol.p, sol.v, sol.f,
                              sol.q̅, sol.p̅, sol.v̅, sol.f̅,
                              int.cache.q̃, int.cache.p̃, int.cache.ṽ, int.cache.f̃,
                              tableau(int).q.c[i], tableau(int).p.c[i])

        for k in eachdim(int)
            int.solver.x[3*(ndims(int)*(i-1)+k-1)+1] = (int.cache.q̃[k] - sol.q[k])/timestep(int)
            int.solver.x[3*(ndims(int)*(i-1)+k-1)+2] = (int.cache.p̃[k] - sol.p[k])/timestep(int)
            int.solver.x[3*(ndims(int)*(i-1)+k-1)+3] = int.cache.ṽ[k]
        end
    end

    for i in 1:pstages(int)
        evaluate!(int.iguess, sol.q, sol.p, sol.v, sol.f,
                              sol.q̅, sol.p̅, sol.v̅, sol.f̅,
                              int.cache.q̃, int.cache.p̃, int.cache.ṽ, int.cache.f̃,
                              tableau(int).q̃.c[i], tableau(int).p̃.c[i])

        for k in eachdim(int)
            int.solver.x[3*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+1] = (int.cache.q̃[k] - sol.q[k])/timestep(int)
            int.solver.x[3*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+2] = (int.cache.p̃[k] - sol.p[k])/timestep(int)
            int.solver.x[3*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+3] = 0
        end
    end

    if isdefined(tableau(int), :λ) && tableau(int).λ.c[1] == 0
        for k in eachdim(int)
            int.solver.x[3*ndims(int)*nstages(int)+3*(k-1)+3] = int.cache.λ[k]
        end
    end

    if isdefined(tableau(int), :d) && length(tableau(int).d) > 0
        for k in eachdim(int)
            int.solver.x[3*ndims(int)*nstages(int)+3*ndims(int)*pstages(int)+k] = 0
        end
    end
end
