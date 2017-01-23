
"Special Partitioned Additive Runge Kutta integrator."
immutable IntegratorSPARK{DT,TT,VT,FT,UT,GT,ΦT} <: Integrator{DT,TT}
    equation::PDAE{DT,TT,VT,FT,UT,GT,ΦT}
    tableau::TableauSPARK{TT}
    Δt::TT

    solver::NonlinearSolver{DT}

    x::Array{DT,1}
    y::Array{DT,1}
    X::Array{DT,2}
    Y::Array{DT,2}
    F::Array{DT,2}
end

"Integrate partitioned DAE with Special Additive Runge Kutta integrator."
function integrate!(int::IntegratorSPARK, s::SolutionPDAE)
    # TODO
end
