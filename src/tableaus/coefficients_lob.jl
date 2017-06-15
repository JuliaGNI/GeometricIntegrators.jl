
function getCoefficientsLobIII2(T=Float64)
    a = Array{T}(
        [[0  0 ]
         [1  0 ]])
    b = Array{T}([1//2, 1//2])
    c = Array{T}([0,    1   ])

    CoefficientsRK(:LobIII, 2, a, b, c)
end

function getCoefficientsLobIII3(T=Float64)
    a = Array{T}([
         [0      0    0 ]
         [1//4  1//4  0 ]
         [0     1     0 ]
        ])
    b = Array{T}([1//6,  2//3,  1//6])
    c = Array{T}([0,     1//2,  1   ])

    CoefficientsRK(:LobIII, 4, a, b, c)
end

function getCoefficientsLobIII4(T=Float64)
    a = Array{T}(@dec128 [
         [      0             0             0     0 ]
         [ (5+√5)/60          1/6   (15-7*√5)/60  0 ]
         [ (5-√5)/60  (15+7*√5)/60          1/6   0 ]
         [      1/6      (5-√5)/12     (5+√5)/12  0 ]
        ])
    b = Array{T}(@dec128 [1/12,      5/12,      5/12,   1/12])
    c = Array{T}(@dec128 [0,    (5-√5)/10, (5+√5)/10,   1   ])

    CoefficientsRK(:LobIII, 6, a, b, c)
end


function getCoefficientsLobIIIA2(T=Float64)
    a = Array{T}([
         [0     0   ]
         [1//2  1//2]
        ])
    b = Array{T}([1//2, 1//2])
    c = Array{T}([0,    1   ])

    CoefficientsRK(:LobIIIA, 2, a, b, c)
end

function getCoefficientsLobIIIA3(T=Float64)
    a = Array{T}([
         [0      0     0    ]
         [5//24  1//3 -1//24]
         [1//6   2//3  1//6 ]
        ])
    b = Array{T}([1//6,  2//3,  1//6])
    c = Array{T}([0,     1//2,  1   ])

    CoefficientsRK(:LobIIIA, 4, a, b, c)
end

function getCoefficientsLobIIIA4(T=Float64)
    a = Array{T}(@dec128 [
         [0             0                0                0           ]
         [(11+√5)/120  (25-   √5)/120  (25-13*√5)/120  (-1+√5)/120]
         [(11-√5)/120  (25+13*√5)/120  (25+   √5)/120  (-1-√5)/120]
         [      1/12            5/12            5/12         1/12 ]
        ])
    b = Array{T}(@dec128 [1/12,      5/12,      5/12,   1/12])
    c = Array{T}(@dec128 [0,    (5-√5)/10, (5+√5)/10,   1    ])

    CoefficientsRK(:LobIIIA, 6, a, b, c)
end


function getCoefficientsLobIIIB2(T=Float64)
    a = Array{T}(
        [[1//2  0]
         [1//2  0]])
    b = Array{T}([1//2, 1//2])
    c = Array{T}([0,    1   ])

    CoefficientsRK(:LobIIIB, 2, a, b, c)
end

function getCoefficientsLobIIIB3(T=Float64)
    a = Array{T}([
         [1//6 -1//6  0   ]
         [1//6  1//3  0   ]
         [1//6  5//6  0   ]
        ])
    b = Array{T}([1//6,  2//3,  1//6])
    c = Array{T}([0,     1//2,  1   ])

    CoefficientsRK(:LobIIIB, 4, a, b, c)
end

function getCoefficientsLobIIIB4(T=Float64)
    a = Array{T}(@dec128 [
         [ 1/12  (-1-   √5)/24   (-1+   √5)/24    0 ]
         [ 1/12  (25+   √5)/120  (25-13*√5)/120   0 ]
         [ 1/12  (25+13*√5)/120  (25-   √5)/120   0 ]
         [ 1/12  (11-   √5)/24   (11+   √5)/24    0 ]
        ])
    b = Array{T}(@dec128 [1/12,       5/12,      5/12,   1/12])
    c = Array{T}(@dec128 [0,     (5-√5)/10, (5+√5)/10,   1    ])

    CoefficientsRK(:LobIIIB, 6, a, b, c)
end


function getCoefficientsLobIIIC2(T=Float64)
    a = Array{T}([
         [1//2 -1//2]
         [1//2  1//2]
        ])
    b = Array{T}([1//2, 1//2])
    c = Array{T}([0,    1   ])

    CoefficientsRK(:LobIIIC, 2, a, b, c)
end

function getCoefficientsLobIIIC3(T=Float64)
    a = Array{T}([
         [1//6 -1//3   1//6 ]
         [1//6  5//12 -1//12]
         [1//6  2//3   1//6 ]
        ])
    b = Array{T}([1//6,  2//3,  1//6])
    c = Array{T}([0,     1//2,  1   ])

    CoefficientsRK(:LobIIIC, 4, a, b, c)
end

function getCoefficientsLobIIIC4(T=Float64)
    a = Array{T}(@dec128 [
         [ 1/12        -√5/12         √5/12   -1/12 ]
         [ 1/12          1/4   (10-7*√5)/60   √5/60 ]
         [ 1/12  (10+7*√5)/60          1/4   -√5/60 ]
         [ 1/12          5/12          5/12    1/12 ]
        ])
    b = Array{T}(@dec128 [1/12,      5/12,      5/12,   1/12])
    c = Array{T}(@dec128 [0,    (5-√5)/10, (5+√5)/10,   1   ])

    CoefficientsRK(:LobIIIC, 6, a, b, c)
end


function getCoefficientsLobIIIF2(T=Float64)
    a = Array{T}(
        [[1//12 -1//12]
         [7//12  5//12]])
    b = Array{T}([1//2, 1//2])
    c = Array{T}([0,    1   ])

    CoefficientsRK(:LobIIIF, 4, a, b, c)
end

function getCoefficientsLobIIIF3(T=Float64)
    a = Array{T}(
        [[1//30  -1//15   1//30]
         [5//24   1//3   -1//24]
         [2//15  11//15   2//15]])
    b = Array{T}([1//6,  2//3,  1//6])
    c = Array{T}([0,     1//2,  1   ])

    CoefficientsRK(:LobIIIF, 6, a, b, c)
end

function getCoefficientsLobIIIF4(T=Float64)
    a = Array{T}(@dec128 [
         [  1/56                 -√5/56           √5/56   -1/56         ]
         [ 37/420+√5/120  5/24-   √5/210  5/24-47*√5/420  -1/210+√5/120 ]
         [ 37/420-√5/120  5/24+47*√5/420  5/24+   √5/210  -1/210-√5/120 ]
         [ 17/168         5/12-   √5/56   5/12+   √5/56   11/168        ]
        ])
    b = Array{T}(@dec128 [1/12,      5/12,      5/12,   1/12])
    c = Array{T}(@dec128 [0,    (5-√5)/10, (5+√5)/10,   1   ])

    CoefficientsRK(:LobIIIF, 8, a, b, c)
end


function getCoefficientsLobIIID2(T=Float64)
    lob  = getCoefficientsLobIII2(T)
    lobC = getCoefficientsLobIIIC2(T)
    CoefficientsRK(:LobIIID, lobC.o, (lob.a + lobC.a)/2, lobC.b, lobC.c)
end

function getCoefficientsLobIIID3(T=Float64)
    lob  = getCoefficientsLobIII3(T)
    lobC = getCoefficientsLobIIIC3(T)
    CoefficientsRK(:LobIIID, lobC.o, (lob.a + lobC.a)/2, lobC.b, lobC.c)
end

function getCoefficientsLobIIID4(T=Float64)
    lob  = getCoefficientsLobIII4(T)
    lobC = getCoefficientsLobIIIC4(T)
    CoefficientsRK(:LobIIID, lobC.o, (lob.a + lobC.a)/2, lobC.b, lobC.c)
end


function getCoefficientsLobIIIE2(T=Float64)
    lobA = getCoefficientsLobIIIA2(T)
    lobB = getCoefficientsLobIIIB2(T)
    CoefficientsRK(:LobIIIE, lobA.o, (lobA.a + lobB.a)/2, lobA.b, lobA.c)
end

function getCoefficientsLobIIIE3(T=Float64)
    lobA = getCoefficientsLobIIIA3(T)
    lobB = getCoefficientsLobIIIB3(T)
    CoefficientsRK(:LobIIIE, lobA.o, (lobA.a + lobB.a)/2, lobA.b, lobA.c)
end

function getCoefficientsLobIIIE4(T=Float64)
    lobA = getCoefficientsLobIIIA4(T)
    lobB = getCoefficientsLobIIIB4(T)
    CoefficientsRK(:LobIIIE, lobA.o, (lobA.a + lobB.a)/2, lobA.b, lobA.c)
end
