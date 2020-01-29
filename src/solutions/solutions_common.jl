
function determine_qdim(equation::Union{SDE,PSDE,SPSDE})
    ns = equation.ns
    ni = equation.ni

    @assert ns ≥ 1
    @assert ni ≥ 1

    if ns==ni==1
        NQ = 2
    elseif ns==1 || ni==1
        NQ = 3
    else
        @error("Both the number of sample paths and the number of initial conditions is larger than one!")
    end

    return NQ
end
