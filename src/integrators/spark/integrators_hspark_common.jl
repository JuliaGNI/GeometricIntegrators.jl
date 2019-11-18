
"Parameters for right-hand side function of Specialised Partitioned Additive Runge-Kutta methods for Variational systems."
mutable struct AbstractParametersHSPARK{IT,DT,TT,D,S,R,P,VT,FT,UT,GT,ϕT,tabType} <: Parameters{DT,TT}
    f_v::VT
    f_f::FT
    f_u::UT
    f_g::GT
    f_ϕ::ϕT

    Δt::TT

    tab::tabType

    t::TT
    q::Vector{DT}
    p::Vector{DT}
    λ::Vector{DT}

    function AbstractParametersHSPARK{IT,DT,D,S,R,P}(f_v::VT, f_f::FT, f_u::UT, f_g::GT, f_ϕ::ϕT, Δt::TT, tableau::tabType) where {IT,DT,TT,D,S,R,P,VT,FT,UT,GT,ϕT,tabType}
        # create solution vectors
        q = zeros(DT,D)
        p = zeros(DT,D)
        λ = zeros(DT,D)

        new{IT,DT,TT,D,S,R,P,VT,FT,UT,GT,ϕT,tabType}(f_v, f_f, f_u, f_g, f_ϕ, Δt, tableau, zero(TT), q, p, λ)
    end
end


function update_params!(params::AbstractParametersHSPARK, cache::IntegratorCacheSPARK)
    # set time for nonlinear solver and copy previous solution
    params.t  = cache.t
    params.q .= cache.q
    params.p .= cache.p
    params.λ .= cache.λ
end


equation(int::AbstractIntegratorHSPARK) = int.equation
timestep(int::AbstractIntegratorHSPARK) = int.params.Δt
tableau(int::AbstractIntegratorHSPARK) = int.tableau
nstages(int::AbstractIntegratorHSPARK) = int.tableau.s
pstages(int::AbstractIntegratorHSPARK) = int.tableau.r


function create_integrator_cache(int::AbstractIntegratorHSPARK{DT,TT}) where {DT,TT}
    IntegratorCacheSPARK{DT, TT, ndims(int), nstages(int), pstages(int)}()
end


function initialize!(int::AbstractIntegratorHSPARK, cache::IntegratorCacheSPARK)
    cache.t̅ = cache.t - timestep(int)

    equation(int).v(cache.t, cache.q, cache.p, cache.v)
    equation(int).f(cache.t, cache.q, cache.p, cache.f)

    initialize!(int.iguess, cache.t, cache.q, cache.p, cache.v, cache.f,
                            cache.t̅, cache.q̅, cache.p̅, cache.v̅, cache.f̅)
end


function initial_guess!(int::AbstractIntegratorHSPARK, cache::IntegratorCacheSPARK)
    for i in 1:nstages(int)
        evaluate!(int.iguess, cache.q, cache.p, cache.v, cache.f,
                              cache.q̅, cache.p̅, cache.v̅, cache.f̅,
                              cache.q̃, cache.p̃, cache.ṽ, cache.f̃,
                              tableau(int).q.c[i], tableau(int).p.c[i])

        for k in 1:ndims(int)
            int.solver.x[2*(ndims(int)*(i-1)+k-1)+1] = (cache.q̃[k] - cache.q[k])/timestep(int)
            int.solver.x[2*(ndims(int)*(i-1)+k-1)+2] = (cache.p̃[k] - cache.p[k])/timestep(int)
        end
    end

    for i in 1:pstages(int)
        evaluate!(int.iguess, cache.q, cache.p, cache.v, cache.f,
                              cache.q̅, cache.p̅, cache.v̅, cache.f̅,
                              cache.q̃, cache.p̃, cache.ṽ, cache.f̃,
                              tableau(int).q̃.c[i], tableau(int).p̃.c[i])

        for k in 1:ndims(int)
            int.solver.x[2*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+1] = (cache.q̃[k] - cache.q[k])/timestep(int)
            int.solver.x[2*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+2] = (cache.p̃[k] - cache.p[k])/timestep(int)
            int.solver.x[2*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+3] = 0
        end
    end

    if tableau(int).λ.c[1] == 0
        for k in 1:ndims(int)
            int.solver.x[2*ndims(int)*nstages(int)+3*(k-1)+3] = cache.λ[k]
        end
    end
end
