
"Special Additive Runge Kutta integrator."
struct IntegratorSARK{DT,TT,VT,UT,ΦT} <: Integrator{DT,TT}
    equation::DAE{DT,TT,VT,UT,ΦT}
    tableau::TableauSARK{TT}

end

"Integrate DAE with Special Additive Runge Kutta integrator."
function integrate!(int::IntegratorSARK, s::SolutionDAE)
    # TODO
end
