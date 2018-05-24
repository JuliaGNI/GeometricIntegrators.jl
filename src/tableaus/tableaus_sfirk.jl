# Tableau for the s-stage Gauss-Lobatto SFIRK method
function getTableauStochasticGLRK(s::Int)

    if s==1
        name = :StochasticGLRK1
    elseif s==2
        name = :StochasticGLRK2
    elseif s==3
        name = :StochasticGLRK3
    elseif s==4
        name = :StochasticGLRK4
    elseif s==5
        name = :StochasticGLRK5
    else s==6
        name = :StochasticGLRK6
    end

    TableauSFIRK(name,getCoefficientsGLRK(s),getCoefficientsGLRK(s))
end
