
function getCoefficientsLobIIIA2()
    a = [[0.0 0.0]
         [0.5 0.5]]
    b = [0.5, 0.5]
    c = [0.0, 1.0]

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

function getCoefficientsLobIIIB2()
    a = [[0.5 0.0]
         [0.5 0.0]]
    b = [0.5, 0.5]
    c = [0.0, 1.0]

    CoefficientsRK(:LobIIIB, 4, a, b, c)
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
