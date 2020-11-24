
import GeometricIntegrators.Quadratures: LobattoLegendreQuadrature, weights
import DecFP: Dec128
import Polynomials


@doc raw"""
The s-stage Lobatto nodes are defined as the roots of the following polynomial of degree $s$:
```math
\frac{d^{s-2}}{dx^{s-2}} \big( (x - x^2)^{s-1} \big) .
```
"""
function get_lobatto_nodes(s)
    if s == 1
        @error "Lobatto nodes for one stage are not defined."
    end
    D(k) = Polynomials.derivative(Polynomials.Polynomial(BigFloat[0, 1, -1])^(k-1), k-2)
    c = sort(real.(Polynomials.roots(D(s))))
    c[begin] = 0; c[end] = 1; c
end

@doc raw"""
The Lobatto weights can be explicitly computed by the formula
```math
b_j = \frac{1}{s (s-1) P_{s-1}(2 c_j - 1)^2} \qquad j = 1 , \, ... , \, s ,
```
where $P_k$ is the $k$th Legendre polynomial, given by
```math
P_k (x) = \frac{1}{k! 2^k} \big( \frac{d^k}{dx^k} (x^2 - 1)^k \big) .
```
"""
function get_lobatto_weights(s)
    if s == 1
        @error "Lobatto weights for one stage are not defined."
    end
    P(k,x) = Polynomials.derivative(Polynomials.Polynomial([-1, 0, 1])^k, k)(x) / factorial(k) / 2^k
    c = get_lobatto_nodes(s)
    b = [ 1 / ( s*(s-1) * P(s-1, 2c[i] - 1)^2 ) for i in 1:s ]
end

@doc raw"""
The Lobatto IIIA coefficients are implicitly given by the so-called simplifying assumption $C(s)$:
```math
\sum \limits_{j=1}^{s} a_{ij} c_{j}^{k-1} = \frac{c_i^k}{k}  \qquad i = 1 , \, ... , \, s , \; k = 1 , \, ... , \, s .
```
"""
function get_lobatto_coefficients_a(s)
    if s == 1
        @error "Lobatto IIIA coefficients for one stage are not defined."
    end

    c = get_lobatto_nodes(s)
    M = [ c[j]^(k-1) for k in 1:s, j in 1:s ]
    
    row(i) = begin
        r = [ c[i]^k / k for k in 1:s ]
        M \ r
    end
    
    vcat([row(i)' for i in 1:s]...)
end

@doc raw"""
The Lobatto IIIB coefficients are implicitly given by the so-called simplifying assumption $D(s)$:
```math
\sum \limits_{i=1}^{s} b_i c_{i}^{k-1} a_{ij} = \frac{b_j}{k} ( 1 - c_j^k)  \qquad j = 1 , \, ... , \, s , \; k = 1 , \, ... , \, s .
```
"""
function get_lobatto_coefficients_b(s)
    if s == 1
        @error "Lobatto IIIB coefficients for one stage are not defined."
    end

    b = get_lobatto_weights(s)
    c = get_lobatto_nodes(s)
    
    M = [ b[i] * c[i]^(k-1) for k in 1:s, i in 1:s ]
    
    row(j) = begin
        r = [ b[j] / k * (1 - c[j]^k) for k in 1:s ]
        M \ r
    end
    
    hcat([row(j) for j in 1:s]...)
end

@doc raw"""
The Lobatto IIIC coefficients are determined by setting $a_{i,1} = b_1$ and
solving the so-called simplifying assumption $C(s-1)$, given by
```math
\sum \limits_{j=1}^{s} a_{ij} c_{j}^{k-1} = \frac{c_i^k}{k}  \qquad i = 1 , \, ... , \, s , \; k = 1 , \, ... , \, s-1 ,
```
for $a_{i,j}$ with $i = 1, ..., s$ and $j = 2, ..., s$.
"""
function get_lobatto_coefficients_c(s)
    if s == 1
        @error "Lobatto IIIC coefficients for one stage are not defined."
    end

    b = get_lobatto_weights(s)
    c = get_lobatto_nodes(s)
    M = [ c[j]^(k-1) for k in 1:s-1, j in 2:s ]
    
    row(i) = begin
        r = [ c[i]^k / k - c[1]^(k-1) * b[1] for k in 1:s-1 ]
        M \ r
    end
    
    hcat(b[1] * ones(s), vcat([row(i)' for i in 1:s]...))
end

@doc raw"""
The Lobatto IIIC̄ coefficients are determined by setting $a_{i,s} = 0$ and
solving the so-called simplifying assumption $C(s-1)$, given by
```math
\sum \limits_{j=1}^{s} a_{ij} c_{j}^{k-1} = \frac{c_i^k}{k}  \qquad i = 1 , \, ... , \, s , \; k = 1 , \, ... , \, s-1 ,
```
for $a_{i,j}$ with $i = 1, ..., s$ and $j = 1, ..., s-1$.
"""
function get_lobatto_coefficients_c̄(s)
    if s == 1
        @error "Lobatto IIIC coefficients for one stage are not defined."
    end

    c = get_lobatto_nodes(s)
    M = [ c[j]^(k-1) for k in 1:s-1, j in 1:s-1 ]
    
    row(i) = begin
        r = [ c[i]^k / k for k in 1:s-1 ]
        M \ r
    end
    
    hcat(vcat([row(i)' for i in 1:s]...), zeros(s))
end


function get_lobatto_coefficients_f(s)
    c = get_lobatto_nodes(s)
    M = [ 1 / (k + j - 1) for k in big.(1:s), j in big.(1:s) ]
    r = [ 1 / s / (s + k) for k in big.(1:s) ]
    α = M \ r
    
    Vₛ = [ c[i]^(j-1) for i in big.(1:s), j in big.(1:s) ]
    Aₛ = zeros(BigFloat, s, s)
    for i in big.(2:s)
       Aₛ[i,i-1] = 1 / (i-1)
    end
    Aₛ[:,s] = α
    
    Vₛ * Aₛ * inv(Vₛ)
end


"Lobatto IIIA coefficients with s stages"
function CoefficientsLobattoIIIA(s, T=Float64)
    CoefficientsRK(T, Symbol("LobattoIIIA($s)"), 2s-2, get_lobatto_coefficients_a(s), get_lobatto_weights(s), get_lobatto_nodes(s))
end

"Lobatto IIIB coefficients with s stages"
function CoefficientsLobattoIIIB(s, T=Float64)
    CoefficientsRK(T, Symbol("LobattoIIIB($s)"), 2s-2, get_lobatto_coefficients_b(s), get_lobatto_weights(s), get_lobatto_nodes(s))
end

"Lobatto IIIC coefficients with s stages"
function CoefficientsLobattoIIIC(s, T=Float64)
    CoefficientsRK(T, Symbol("LobattoIIIC($s)"), 2s-2, get_lobatto_coefficients_c(s), get_lobatto_weights(s), get_lobatto_nodes(s))
end

"Lobatto IIIC̄ coefficients with s stages"
function CoefficientsLobattoIIIC̄(s, T=Float64)
    CoefficientsRK(T, Symbol("LobattoIIIC̄($s)"), 2s-2, get_lobatto_coefficients_c̄(s), get_lobatto_weights(s), get_lobatto_nodes(s))
end

"Lobatto IIID coefficients with s stages"
function CoefficientsLobattoIIID(s, T=Float64)
    CoefficientsRK(T, Symbol("LobattoIIID($s)"), 2s-2, (get_lobatto_coefficients_c(s) .+ get_lobatto_coefficients_c̄(s)) ./ 2, get_lobatto_weights(s), get_lobatto_nodes(s))
end

"Lobatto IIIE coefficients with s stages"
function CoefficientsLobattoIIIE(s, T=Float64)
    CoefficientsRK(T, Symbol("LobattoIIIE($s)"), 2s-2, (get_lobatto_coefficients_a(s) .+ get_lobatto_coefficients_b(s)) ./ 2, get_lobatto_weights(s), get_lobatto_nodes(s))
end

"Lobatto IIIF coefficients with s stages"
function CoefficientsLobattoIIIF(s, T=Float64)
    CoefficientsRK(T, Symbol("LobattoIIIF($s)"), 2s, get_lobatto_coefficients_f(s), get_lobatto_weights(s), get_lobatto_nodes(s))
end

"Lobatto IIIG coefficients with s stages"
function CoefficientsLobattoIIIG(s, T=Float64)
    symplecticize(CoefficientsLobattoIIIF(s, BigFloat); name=Symbol("LobattoIIIG($s)"), T=T)
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
        @error("Number of stages s=$(s) and σ=$(σ) is not supported.")
    end

    CoefficientsIRK{T}(:LobattoIIIIS, s^2, s, σ, a, get_lobatto_weights(σ), get_lobatto_nodes(σ))
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
        @error("We don't have a d vector for s=$(s) stages.")
    end
    return d
end

function get_lobatto_ω_matrix(s)
    as = CoefficientsLobattoIIIA(s).a[2:s,1:s]
    es = zeros(s)
    es[s] = 1

    Q = vcat( hcat(as, zeros(s-1)), hcat(zeros(s)', 1) )
    L = vcat( as, es' )
    ω = inv(L) * Q
    
    return ω
end
