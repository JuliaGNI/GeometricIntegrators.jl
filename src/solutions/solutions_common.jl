
function determine_qdim(equation::Union{SDE,PSDE,SPSDE})
    nd = equation.d
    ns = equation.ns
    ni = equation.ni

    if nd==ns==ni==1
        NQ = 1
    elseif ns==ni==1
        NQ = 2
    elseif ns==1 || ni==1
        NQ = 3
    else
        NQ = 4
    end

    return NQ
end
