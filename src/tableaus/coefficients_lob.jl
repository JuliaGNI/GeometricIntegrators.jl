

function get_lobatto_nodes(s)
    if s == 2
        c = @dec128 [0, 1]
    elseif s == 3
        c = @dec128 [0,    1/2,    1]
    elseif s == 4
        c = @dec128 [0,    (5-√5)/10,    (5+√5)/10,    1]
    elseif s == 5
        c = @dec128 [0,    1/2-√21/14,    1/2,    1/2+√21/14,    1]
    else
        @error "Lobatto nodes for " * string(s) * " stages not implemented."
    end

    return c
end

function get_lobatto_weights(s)
    if s == 2
        b = @dec128 [1/2, 1/2]
    elseif s == 3
        b = @dec128 [1/6,  2/3,  1/6]
    elseif s == 4
        b = @dec128 [1/12,   5/12,   5/12,   1/12]
    elseif s == 5
        b = @dec128 [1/20,   49/180,   16/45,   49/180,   1/20]
    else
        @error "Lobatto weights for " * string(s) * " stages not implemented."
    end

    return b
end


function getCoefficientsLobIII2(T=Float64)
    a = @dec128 [
            [0  0]
            [1  0]
        ]

    CoefficientsRK(T, :LobIII2, 2, a, get_lobatto_weights(2), get_lobatto_nodes(2))
end

function getCoefficientsLobIII3(T=Float64)
    a = @dec128 [
            [0    0    0 ]
            [1/4  1/4  0 ]
            [0    1    0 ]
        ]

    CoefficientsRK(T, :LobIII3, 4, a, get_lobatto_weights(3), get_lobatto_nodes(3))
end

function getCoefficientsLobIII4(T=Float64)
    a = @dec128 [
         [      0             0             0     0 ]
         [ (5+√5)/60          1/6   (15-7*√5)/60  0 ]
         [ (5-√5)/60  (15+7*√5)/60          1/6   0 ]
         [      1/6      (5-√5)/12     (5+√5)/12  0 ]
        ]

    CoefficientsRK(T, :LobIII4, 6, a, get_lobatto_weights(4), get_lobatto_nodes(4))
end


function getCoefficientsLobIII5(T=Float64)
    a = @dec128 [
         [ 0     0                0              0                0 ]
         [ 1/14            1/9    (13-3*√21)/63  (14- 3*√21)/126  0 ]
         [ 1/32  (91+21*√21)/576          11/72  (91-21*√21)/576  0 ]
         [ 1/14  (14+ 3*√21)/126  (13+3*√21)/63            1/9    0 ]
         [ 0               7/18            2/9             7/18   0 ]
        ]

    CoefficientsRK(T, :LobIII5, 8, a, get_lobatto_weights(5), get_lobatto_nodes(5))
end


function getCoefficientsLobIIIA2(T=Float64)
    a = @dec128 [
            [0     0   ]
            [1/2   1/2 ]
        ]

    CoefficientsRK(T, :LobIIIA2, 2, a, get_lobatto_weights(2), get_lobatto_nodes(2))
end

function getCoefficientsLobIIIA3(T=Float64)
    a = @dec128 [
            [0      0     0    ]
            [5/24   1/3  -1/24 ]
            [1/6    2/3   1/6  ]
        ]

    CoefficientsRK(T, :LobIIIA3, 4, a, get_lobatto_weights(3), get_lobatto_nodes(3))
end

function getCoefficientsLobIIIA4(T=Float64)
    a = @dec128 [
         [0            0               0               0          ]
         [(11+√5)/120  (25-   √5)/120  (25-13*√5)/120  (-1+√5)/120]
         [(11-√5)/120  (25+13*√5)/120  (25+   √5)/120  (-1-√5)/120]
         [      1/12            5/12            5/12         1/12 ]
        ]

    CoefficientsRK(T, :LobIIIA4, 6, a, get_lobatto_weights(4), get_lobatto_nodes(4))
end


function getCoefficientsLobIIIA5(T=Float64)
    a = @dec128 [
         [0                 0                   0                  0                   0               ]
         [(119+3*√21)/1960  (343-  9*√21)/2520  (392-96*√21)/2205  (343- 69*√21)/2520  (-21+3*√21)/1960]
         [         13/320   (392+105*√21)/2880             8/45    (392-105*√21)/2880            3/320 ]
         [(119-3*√21)/1960  (343+ 69*√21)/2520  (392+96*√21)/2205  (343+  9*√21)/2520  (-21-3*√21)/1960]
         [          1/20               49/180             16/45               49/180             1/20  ]
        ]

    CoefficientsRK(T, :LobIIIA5, 8, a, get_lobatto_weights(5), get_lobatto_nodes(5))
end


function getCoefficientsLobIIIB2(T=Float64)
    a = @dec128 [
            [1/2  0]
            [1/2  0]
        ]

    CoefficientsRK(T, :LobIIIB2, 2, a, get_lobatto_weights(2), get_lobatto_nodes(2))
end

function getCoefficientsLobIIIB3(T=Float64)
    a = @dec128 [
            [1/6  -1/6   0   ]
            [1/6   1/3   0   ]
            [1/6   5/6   0   ]
        ]

    CoefficientsRK(T, :LobIIIB3, 4, a, get_lobatto_weights(3), get_lobatto_nodes(3))
end

function getCoefficientsLobIIIB4(T=Float64)
    a = @dec128 [
         [ 1/12  (-1-   √5)/24   (-1+   √5)/24    0 ]
         [ 1/12  (25+   √5)/120  (25-13*√5)/120   0 ]
         [ 1/12  (25+13*√5)/120  (25-   √5)/120   0 ]
         [ 1/12  (11-   √5)/24   (11+   √5)/24    0 ]
        ]

    CoefficientsRK(T, :LobIIIB4, 6, a, get_lobatto_weights(4), get_lobatto_nodes(4))
end


function getCoefficientsLobIIIB5(T=Float64)
    a = @dec128 [
         [ 1/20  ( -7-   √21)/120             1/15   ( -7+   √21)/120    0 ]
         [ 1/20  (343+ 9*√21)/2520  (56-15*√21)/315  (343-69*√21)/2520   0 ]
         [ 1/20  ( 49+12*√21)/360             8/45   ( 49-12*√21)/360    0 ]
         [ 1/20  (343+69*√21)/2520  (56+15*√21)/315  (343- 9*√21)/2520   0 ]
         [ 1/20  (119- 3*√21)/360            13/45   (119+ 3*√21)/360    0 ]
        ]

    CoefficientsRK(T, :LobIIIB5, 8, a, get_lobatto_weights(5), get_lobatto_nodes(5))
end


function getCoefficientsLobIIIC2(T=Float64)
    a = @dec128 [
            [1/2  -1/2 ]
            [1/2   1/2 ]
        ]

    CoefficientsRK(T, :LobIIIC2, 2, a, get_lobatto_weights(2), get_lobatto_nodes(2))
end

function getCoefficientsLobIIIC3(T=Float64)
    a = @dec128 [
            [1/6  -1/3    1/6  ]
            [1/6   5/12  -1/12 ]
            [1/6   2/3    1/6  ]
        ]

    CoefficientsRK(T, :LobIIIC3, 4, a, get_lobatto_weights(3), get_lobatto_nodes(3))
end

function getCoefficientsLobIIIC4(T=Float64)
    a = @dec128 [
         [ 1/12        -√5/12         √5/12   -1/12 ]
         [ 1/12          1/4   (10-7*√5)/60   √5/60 ]
         [ 1/12  (10+7*√5)/60          1/4   -√5/60 ]
         [ 1/12          5/12          5/12    1/12 ]
        ]

    CoefficientsRK(T, :LobIIIC4, 6, a, get_lobatto_weights(4), get_lobatto_nodes(4))
end


function getCoefficientsLobIIIC5(T=Float64)
    a = @dec128 [
         [ 1/20            -21/180             2/15             -21/180     1/20  ]
         [ 1/20             29/180   (47-15*√21)/315  (203- 30*√21)/1260   -3/140 ]
         [ 1/20  (329+105*√21)/2880           73/360  (329-105*√21)/2880    3/160 ]
         [ 1/20  (203+ 30*√21)/1260  (47+15*√21)/315             29/180    -3/140 ]
         [ 1/20             49/180            16/45              49/180     1/20  ]
        ]

    CoefficientsRK(T, :LobIIIC5, 8, a, get_lobatto_weights(5), get_lobatto_nodes(5))
end


function getCoefficientsLobIIIF2(T=Float64)
    a = @dec128 [
            [1/12  -1/12 ]
            [7/12   5/12 ]
        ]

    CoefficientsRK(T, :LobIIIF2, 4, a, get_lobatto_weights(2), get_lobatto_nodes(2))
end

function getCoefficientsLobIIIF3(T=Float64)
    a = @dec128 [
            [1/30   -1/15    1/30 ]
            [5/24    1/3    -1/24 ]
            [2/15   11/15    2/15 ]
        ]

    CoefficientsRK(T, :LobIIIF3, 6, a, get_lobatto_weights(3), get_lobatto_nodes(3))
end

function getCoefficientsLobIIIF4(T=Float64)
    a = @dec128 [
            [  1/56                 -√5/56           √5/56   -1/56         ]
            [ 37/420+√5/120  5/24-   √5/210  5/24-47*√5/420  -1/210+√5/120 ]
            [ 37/420-√5/120  5/24+47*√5/420  5/24+   √5/210  -1/210-√5/120 ]
            [ 17/168         5/12-   √5/56   5/12+   √5/56   11/168        ]
        ]

    CoefficientsRK(T, :LobIIIF4, 8, a, get_lobatto_weights(4), get_lobatto_nodes(4))
end


function getCoefficientsLobIIID2(T=Float64)
    lob  = getCoefficientsLobIII2(Dec128)
    lobC = getCoefficientsLobIIIC2(Dec128)
    CoefficientsRK(T, :LobIIID2, lobC.o, (lob.a + lobC.a)/2, lobC.b, lobC.c)
end

function getCoefficientsLobIIID3(T=Float64)
    lob  = getCoefficientsLobIII3(Dec128)
    lobC = getCoefficientsLobIIIC3(Dec128)
    CoefficientsRK(T, :LobIIID3, lobC.o, (lob.a + lobC.a)/2, lobC.b, lobC.c)
end

function getCoefficientsLobIIID4(T=Float64)
    lob  = getCoefficientsLobIII4(Dec128)
    lobC = getCoefficientsLobIIIC4(Dec128)
    CoefficientsRK(T, :LobIIID4, lobC.o, (lob.a + lobC.a)/2, lobC.b, lobC.c)
end

function getCoefficientsLobIIID5(T=Float64)
    lob  = getCoefficientsLobIII5(Dec128)
    lobC = getCoefficientsLobIIIC5(Dec128)
    CoefficientsRK(T, :LobIIID5, lobC.o, (lob.a + lobC.a)/2, lobC.b, lobC.c)
end


function getCoefficientsLobIIIE2(T=Float64)
    lobA = getCoefficientsLobIIIA2(Dec128)
    lobB = getCoefficientsLobIIIB2(Dec128)
    CoefficientsRK(T, :LobIIIE2, lobA.o, (lobA.a + lobB.a)/2, lobA.b, lobA.c)
end

function getCoefficientsLobIIIE3(T=Float64)
    lobA = getCoefficientsLobIIIA3(Dec128)
    lobB = getCoefficientsLobIIIB3(Dec128)
    CoefficientsRK(T, :LobIIIE3, lobA.o, (lobA.a + lobB.a)/2, lobA.b, lobA.c)
end

function getCoefficientsLobIIIE4(T=Float64)
    lobA = getCoefficientsLobIIIA4(Dec128)
    lobB = getCoefficientsLobIIIB4(Dec128)
    CoefficientsRK(T, :LobIIIE4, lobA.o, (lobA.a + lobB.a)/2, lobA.b, lobA.c)
end

function getCoefficientsLobIIIE5(T=Float64)
    lobA = getCoefficientsLobIIIA5(Dec128)
    lobB = getCoefficientsLobIIIB5(Dec128)
    CoefficientsRK(T, :LobIIIE5, lobA.o, (lobA.a + lobB.a)/2, lobA.b, lobA.c)
end


function getCoefficientsLobIIIG2(T=Float64)
    symplecticize(getCoefficientsLobIIIF2(Dec128); name=:LobIIIG2, T=T)
end

function getCoefficientsLobIIIG3(T=Float64)
    symplecticize(getCoefficientsLobIIIF3(Dec128); name=:LobIIIG3, T=T)
end

function getCoefficientsLobIIIG4(T=Float64)
    symplecticize(getCoefficientsLobIIIF4(Dec128); name=:LobIIIG4, T=T)
end


function getCoefficientsLobIII(s, T=Float64)
    if s == 2
        getCoefficientsLobIII2(T)
    elseif s == 3
        getCoefficientsLobIII3(T)
    elseif s == 4
        getCoefficientsLobIII4(T)
    elseif s == 5
        getCoefficientsLobIII5(T)
    else
        @error "Lobatto III Tableau with " * string(s) * " stages not implemented."
    end
end

function getCoefficientsLobIIIA(s, T=Float64)
    if s == 2
        getCoefficientsLobIIIA2(T)
    elseif s == 3
        getCoefficientsLobIIIA3(T)
    elseif s == 4
        getCoefficientsLobIIIA4(T)
    elseif s == 5
        getCoefficientsLobIIIA5(T)
    else
        @error "Lobatto IIIA Tableau with " * string(s) * " stages not implemented."
    end
end

function getCoefficientsLobIIIB(s, T=Float64)
    if s == 2
        getCoefficientsLobIIIB2(T)
    elseif s == 3
        getCoefficientsLobIIIB3(T)
    elseif s == 4
        getCoefficientsLobIIIB4(T)
    elseif s == 5
        getCoefficientsLobIIIB5(T)
    else
        @error "Lobatto IIIB Tableau with " * string(s) * " stages not implemented."
    end
end

function getCoefficientsLobIIIC(s, T=Float64)
    if s == 2
        getCoefficientsLobIIIC2(T)
    elseif s == 3
        getCoefficientsLobIIIC3(T)
    elseif s == 4
        getCoefficientsLobIIIC4(T)
    elseif s == 5
        getCoefficientsLobIIIC5(T)
    else
        @error "Lobatto IIIC Tableau with " * string(s) * " stages not implemented."
    end
end

function getCoefficientsLobIIID(s, T=Float64)
    if s == 2
        getCoefficientsLobIIID2(T)
    elseif s == 3
        getCoefficientsLobIIID3(T)
    elseif s == 4
        getCoefficientsLobIIID4(T)
    elseif s == 5
        getCoefficientsLobIIID5(T)
    else
        @error "Lobatto IIID Tableau with " * string(s) * " stages not implemented."
    end
end

function getCoefficientsLobIIIE(s, T=Float64)
    if s == 2
        getCoefficientsLobIIIE2(T)
    elseif s == 3
        getCoefficientsLobIIIE3(T)
    elseif s == 4
        getCoefficientsLobIIIE4(T)
    elseif s == 5
        getCoefficientsLobIIIE5(T)
    else
        @error "Lobatto IIIE Tableau with " * string(s) * " stages not implemented."
    end
end

function getCoefficientsLobIIIF(s, T=Float64)
    if s == 2
        getCoefficientsLobIIIF2(T)
    elseif s == 3
        getCoefficientsLobIIIF3(T)
    elseif s == 4
        getCoefficientsLobIIIF4(T)
    else
        @error "Lobatto IIIF Tableau with " * string(s) * " stages not implemented."
    end
end

function getCoefficientsLobIIIG(s, T=Float64)
    if s == 2
        getCoefficientsLobIIIG2(T)
    elseif s == 3
        getCoefficientsLobIIIG3(T)
    elseif s == 4
        getCoefficientsLobIIIG4(T)
    else
        @error "Lobatto IIIG Tableau with " * string(s) * " stages not implemented."
    end
end


function get_lobatto_projective_stage(s, T=Float64)
    if s == 1
        a = reshape(T[0,1], (2,1))
    elseif s == 2
        a = T[ 0     0
               1//2  1//2 ]
    elseif s == 3
        a = T[ 0      0      0
               5//18  8//18  5//18 ]
    elseif s == 4
        a = T[ 0            0            0            0
               1//4-√30/72  1//4+√30/72  1//4+√30/72  1//4-√30/72 ]
    elseif s == 5
        a = T[ 0            0            0            0            0
               (-5797*√70 - 19635*√6 + 3297*√105 + 34069)/(900*(-31*√70 - 105*√6 + 12*√105 + 124))  3*(-155*√70 - 525*√6 + 21*√105 + 217)/(100*(-31*√70 - 105*√6 + 12*√105 + 124))  64/225  3*(-155*√70 - 525*√6 + 21*√105 + 217)/(100*(-31*√70 - 105*√6 + 12*√105 + 124))  (-5797*√70 - 19635*√6 + 3297*√105 + 34069)/(900*(-31*√70 - 105*√6 + 12*√105 + 124)) ]
            #    0.1184634425280945437571320203599586813216300011062070077914139441108586442015225  0.2393143352496832340206457574178190964561477766715707699863638336669191335762568  64//225  0.2393143352496832340206457574178190964561477766715707699863638336669191335762568  0.1184634425280945437571320203599586813216300011062070077914139441108586442015225
    else
        @error "Lobatto projective coefficients for " * string(s) * " stages not implemented."
    end

    return a
end


function get_lobatto_interstage_coefficients(s, σ=s+1, T=Float64)
    if s == 1 && σ == 2
        a = reshape(Array{Dec128}(@dec128 [
                [0]
                [1]
            ]), σ, s)
    elseif s == 2 && σ == 3
        a = @dec128 [
                [0         0       ]
                [1/4+√3/8  1/4-√3/8]
                [1/2       1/2     ]
            ]
    elseif s == 3 && σ == 4
        a = @dec128 [
                [0                     0              0                 ]
                [5/36-√5/180+√15/30    2/9-4*√5/45    5/36-√15/30-√5/180]
                [5/36+√5/180+√15/30    2/9+4*√5/45    5/36+√5/180-√15/30]
                [5/18                  4/9            5/18              ]
            ]
    elseif s == 4 && σ == 5
        a = @dec128[
                [ 0  0  0  0 ]
                [ 0.1591671294900345358547852036503968244037099929876660089701433805420016628303049  0.01565551541810189985834126767767053597372403599533987485875443079743596295789264  -0.002649931830622319490548503812025451589560291885905677307276227674989121280824921  0.0005004515684973118782758043605289134275222217051875147207182646279180365249293946 ]
                [ 0.176105937938662021225385981377510389121610166295011341406551070207182545793024   0.3049159394838621526624258370570061164298248408582482915516281042450345033436905    0.02115663794741091865104218833199417995250081122480493820210723599558956292335403  -0.002178515369935092538854006766510685503935818378064571160286410447806612060067542  ]
                [ 0.1734269710002296168082561702504707901901521262117592555255463951314578972080256  0.3287225092618953908040165292010257479718859439689589070610115679156131875478712    0.3104170620131711714551267577113297604086016160877133548949809094431881033091524    0.01476029307869239283174677096060287921396435492928076127612127921737427090265032   ]
                [ 0.1739274225687269286865319746109997036176743479169467702462646597593759337329541  0.3260725774312730713134680253890002963823256520830532297537353402406240662670459    0.3260725774312730713134680253890002963823256520830532297537353402406240662670459    0.1739274225687269286865319746109997036176743479169467702462646597593759337329541    ]
            ]            
    else
        @error("Number of stages s=" * string(s) * " and σ=" * string(σ) * " is not supported.")
    end

    CoefficientsIRK{T}(:LobIIIIS, s^2, s, σ, a, get_lobatto_weights(σ), get_lobatto_nodes(σ))
end


function get_lobatto_d_vector(s)
    if s == 2
        d = [+1.0, -1.0]
    elseif s == 3
        d = [+1.0, -2.0, +1.0]
    elseif s == 4
        d = [+1.0, -√5, +√5, -1.0]
    elseif s == 5
        d = [+3.0, -7.0, +8.0, -7.0, +3.0]
    else
        @error("We don't have a d vector for s=" * string(s) * " stages.")
    end
    return d
end

function get_lobatto_ω_matrix(s)
    as = getCoefficientsLobIIIA(s).a[2:s,1:s]
    es = zeros(s)
    es[s] = 1

    Q = vcat( hcat(as, zeros(s-1)), hcat(zeros(s)', 1) )
    L = vcat( as, es' )
    ω = inv(L) * Q
    
    return ω
end
