
function getCoefficientsLobIII2()
    a = Array{Float64}(
        [[0  0 ]
         [1  0 ]])
    b = Array{Float64}([1//2, 1//2])
    c = Array{Float64}([0,    1   ])

    CoefficientsRK(:LobIII, 2, a, b, c)
end

function getCoefficientsLobIII3()
    a = Array{Float64}(
        [[0      0    0 ]
         [1//4  1//4  0 ]
         [0     1     0 ]])
    b = Array{Float64}([1//6,  2//3,  1//6])
    c = Array{Float64}([0,     1//2,  1   ])

    CoefficientsRK(:LobIII, 4, a, b, c)
end

function getCoefficientsLobIII4()
    a = Array{Float64}(
        [[       0             0             0    0 ]
         [ (5+√5)/60          1/6   (15-7*√5)/60  0 ]
         [ (5-√5)/60  (15+7*√5)/60          1/6   0 ]
         [      1/6      (5-√5)/12     (5+√5)/12  0 ]])
    b = Array{Float64}([1/12,      5/12,      5/12,   1/12])
    c = Array{Float64}([0,    (5-√5)/10, (5+√5)/10,   1    ])

    CoefficientsRK(:LobIII, 6, a, b, c)
end


function getCoefficientsLobIIIA2()
    a = Array{Float64}(
        [[0     0   ]
         [1//2  1//2]])
    b = Array{Float64}([1//2, 1//2])
    c = Array{Float64}([0,    1   ])

    CoefficientsRK(:LobIIIA, 2, a, b, c)
end

function getCoefficientsLobIIIA3()
    a = Array{Float64}(
        [[0      0     0    ]
         [5//24  1//3 -1//24]
         [1//6   2//3  1//6 ]])
    b = Array{Float64}([1//6,  2//3,  1//6])
    c = Array{Float64}([0,     1//2,  1   ])

    CoefficientsRK(:LobIIIA, 4, a, b, c)
end

function getCoefficientsLobIIIA4()
    a = Array{Float64}(
        [[0             0                0                0           ]
         [(11+√5)/120  (25-   √5)/120  (25-13*√5)/120  (-1+√5)/120]
         [(11-√5)/120  (25+13*√5)/120  (25+   √5)/120  (-1-√5)/120]
         [      1/12            5/12            5/12         1/12 ]])
    b = Array{Float64}([1/12,      5/12,      5/12,   1/12])
    c = Array{Float64}([0,    (5-√5)/10, (5+√5)/10,   1    ])

    CoefficientsRK(:LobIIIA, 6, a, b, c)
end


function getCoefficientsLobIIIB2()
    a = Array{Float64}(
        [[1//2  0]
         [1//2  0]])
    b = Array{Float64}([1//2, 1//2])
    c = Array{Float64}([0,    1   ])

    CoefficientsRK(:LobIIIB, 2, a, b, c)
end

function getCoefficientsLobIIIB3()
    a = Array{Float64}(
        [[1//6 -1//6  0   ]
         [1//6  1//3  0   ]
         [1//6  5//6  0   ]])
    b = Array{Float64}([1//6,  2//3,  1//6])
    c = Array{Float64}([0,     1//2,  1   ])

    CoefficientsRK(:LobIIIB, 4, a, b, c)
end

function getCoefficientsLobIIIB4()
    a = Array{Float64}(
        [[ 1/12  (-1-   √5)/24   (-1+   √5)/24    0 ]
         [ 1/12  (25+   √5)/120  (25-13*√5)/120   0 ]
         [ 1/12  (25+13*√5)/120  (25-   √5)/120   0 ]
         [ 1/12  (11-   √5)/24   (11+   √5)/24    0 ]])
    b = Array{Float64}([1/12,       5/12,      5/12,   1/12])
    c = Array{Float64}([0,     (5-√5)/10, (5+√5)/10,   1   ])

    CoefficientsRK(:LobIIIB, 6, a, b, c)
end


function getCoefficientsLobIIIC2()
    a = Array{Float64}(
        [[1//2 -1//2]
         [1//2  1//2]])
    b = Array{Float64}([1//2, 1//2])
    c = Array{Float64}([0,    1   ])

    CoefficientsRK(:LobIIIC, 2, a, b, c)
end

function getCoefficientsLobIIIC3()
    a = Array{Float64}(
        [[1//6 -1//3   1//6 ]
         [1//6  5//12 -1//12]
         [1//6  2//3   1//6 ]])
    b = Array{Float64}([1//6,  2//3,  1//6])
    c = Array{Float64}([0,     1//2,  1   ])

    CoefficientsRK(:LobIIIC, 4, a, b, c)
end

function getCoefficientsLobIIIC4()
    a = Array{Float64}(
        [[ 1/12        -√5/12         √5/12   -1/12 ]
         [ 1/12          1/4   (10-7*√5)/60   √5/60 ]
         [ 1/12  (10+7*√5)/60          1/4   -√5/60 ]
         [ 1/12          5/12          5/12    1/12 ]])
    b = Array{Float64}([1//12,      5/12,      5/12,   1/12])
    c = Array{Float64}([0,     (5-√5)/10, (5+√5)/10,   1   ])

    CoefficientsRK(:LobIIIC, 6, a, b, c)
end


function getCoefficientsLobIIID2()
    lob  = getCoefficientsLobIII2()
    lobC = getCoefficientsLobIIIC2()
    CoefficientsRK(:LobIIID, lobC.o, 0.5*(lob.a + lobC.a), lobC.b, lobC.c)
end

function getCoefficientsLobIIID3()
    lob  = getCoefficientsLobIII3()
    lobC = getCoefficientsLobIIIC3()
    CoefficientsRK(:LobIIID, lobC.o, 0.5*(lob.a + lobC.a), lobC.b, lobC.c)
end

function getCoefficientsLobIIID4()
    lob  = getCoefficientsLobIII4()
    lobC = getCoefficientsLobIIIC4()
    CoefficientsRK(:LobIIID, lobC.o, 0.5*(lob.a + lobC.a), lobC.b, lobC.c)
end


function getCoefficientsLobIIIE2()
    lobA = getCoefficientsLobIIIA2()
    lobB = getCoefficientsLobIIIB2()
    CoefficientsRK(:LobIIIE, lobA.o, 0.5*(lobA.a + lobB.a), lobA.b, lobA.c)
end

function getCoefficientsLobIIIE3()
    lobA = getCoefficientsLobIIIA3()
    lobB = getCoefficientsLobIIIB3()
    CoefficientsRK(:LobIIIE, lobA.o, 0.5*(lobA.a + lobB.a), lobA.b, lobA.c)
end

function getCoefficientsLobIIIE4()
    lobA = getCoefficientsLobIIIA4()
    lobB = getCoefficientsLobIIIB4()
    CoefficientsRK(:LobIIIE, lobA.o, 0.5*(lobA.a + lobB.a), lobA.b, lobA.c)
end
