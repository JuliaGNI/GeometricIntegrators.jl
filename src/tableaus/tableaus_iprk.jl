
"Gauss-Legendre Runge-Kutta"
function getTableauIPGLRK(s::Int)
    TableauIPRK(Symbol("IPGLRK", s), 2^s, getCoefficientsGLRK(s))
end
