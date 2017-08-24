
"Discontinuous Galerkin Variational Integrator."
struct IntegratorDGVI{DT,TT,ΘT,FT,GT,VT,FPT,ST,IT,BT<:Basis,D,S,R} <: Integrator{DT,TT}
    equation::IODE{DT,TT,ΘT,FT,GT,VT}

    basis::BT
    quadrature::Quadrature{TT,R}

    Δt::TT

    params::FPT
    solver::ST
    iguess::InitialGuessPODE{DT,TT,VT,FT,IT}

    q::Vector{DT}
    p::Vector{DT}

    # cache::NonlinearFunctionCacheDGVI{DT}
end
