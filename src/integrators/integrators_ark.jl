
"Additive Runge Kutta integrator."
immutable IntegratorARK{DT,TT,VT,UT,ΦT} <: Integrator{DT,TT}
    equation::DAE{DT,TT,VT,UT,ΦT}
    tableau::TableauARK{TT}

end

"Integrate DAE with Additive Runge Kutta integrator."
function integrate!(int::IntegratorARK, s::SolutionDAE)
    # TODO
end
