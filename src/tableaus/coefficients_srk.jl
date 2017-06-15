
function getCoefficientsSRK3(T=Float64)
    a = Array{T}(@dec128 [
         [5/36         2/9        5/36-√15/10]
         [5/36         2/9        5/36       ]
         [5/36+√15/10  2/9        5/36       ]
        ])
    b = Array{T}(@dec128 [5/18,        4/9,       5/18       ])
    c = Array{T}(@dec128 [1/2-√15/10,  1/2,       1/2+√15/10 ])
    o = 4

    CoefficientsRK(:srk3, o, a, b, c)
end
