

function getTableauLieA(T=Float64)
    a = Array{T}([ 1 ])
    b = Array{T}([ 0 ])
    TableauSplittingNS(:LieTrotterSplittingA, 1, a, b)
end


function getTableauLieB(T=Float64)
    a = Array{T}([ 0 ])
    b = Array{T}([ 1 ])
    TableauSplittingNS(:LieTrotterSplittingB, 1, a, b)
end


function getTableauStrang(T=Float64)
    a = Array{T}([ 1//2 ])
    b = Array{T}([ 1//2 ])
    TableauSplittingNS(:StrangSplitting, 2, a, b)
end


function getTableauMcLachlan2(T=Float64; α=0.1932)
    a = Array{T}([ α, 0.5 - α ])
    b = Array{T}([ 0.5 - α, α ])
    TableauSplittingNS(:McLachlanSplitting, 2, a, b)
end


function getTableauMcLachlan4(T=Float64)
    a = Array{T}(@dec128 [ (146 +  5*√19) / 540,
                           ( -2 + 10*√19) / 135,
                           1/5,
                           (-23 - 20*√19) / 270,
                           ( 14 -    √19) / 108])
    TableauSplittingNS(:McLachlanSplitting, 4, a, a[end:-1:1])
end


function getTableauTripleJump(T=Float64)
    fac = @dec128 2^(1/3)
    den = @dec128 1/(2-fac)
    a = Array{T}([ den, -fac*den ])
    TableauSplittingSS(:TripleJumpSplitting, 4, a)
end


function getTableauSuzukiFractal(T=Float64)
    fac = @dec128 4^(1/3)
    den = @dec128 1/(4-fac)
    a = Array{T}([ den, den, -fac*den ])
    TableauSplittingSS(:SuzukiFractalSplitting, 4, a)
end
