
function getCoefficientsSRK3()
    a = [[5/36         2/9        5/36-√15/10]
         [5/36         2/9        5/36       ]
         [5/36+√15/10  2/9        5/36       ]]
    b =  [5/18,        4/9,       5/18       ]
    c =  [1/2-√15/10,  1/2,       1/2+√15/10 ]
    o = 4

    CoefficientsRK(:srk3, o, a, b, c)
end
