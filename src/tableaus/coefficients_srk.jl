
function CoefficientsSRK3(T=Float64)
    a = @dec128 [
         [5/36         2/9        5/36-√15/10]
         [5/36         2/9        5/36       ]
         [5/36+√15/10  2/9        5/36       ]
        ]
    b = @dec128 [5/18,        4/9,       5/18       ]
    c = @dec128 [1/2-√15/10,  1/2,       1/2+√15/10 ]
    o = 4

    CoefficientsRK(T, :SRK3, o, a, b, c)
end
