
function getCoefficientsRadIIA2(T=Float64)
    a = [[5//12  -1//12]
         [3//4    1//4 ]]
    b = [3//4, 1//4]
    c = [1//3, 1//1]

    CoefficientsRK(T, :RadIIA2, 3, a, b, c)
end

function getCoefficientsRadIIA3(T=Float64)
    a = @dec128 [
          [11/45 -  7*√6/360    37/225-169*√6/1800   -2/225+√6/75]
          [37/225+169*√6/1800   11/45 +  7*√6/360    -2/225-√6/75]
          [ 4/9  -    √6/36      4/9  +    √6/36      1/9        ]
        ]
    b = @dec128 [4/9-√6/36,   4/9+√6/36,  1/9 ]
    c = @dec128 [2/5-√6/10,   2/5+√6/10,  1   ]

    CoefficientsRK(T, :RadIIA4, 5, a, b, c)
end
