
import LinearAlgebra
import Polynomials

import GeometricIntegrators.Quadratures: GaussLegendreQuadrature


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


@doc raw"""
The projective Lobatto-GLRK coefficients are implicitly given by
```math
\sum \limits_{j=1}^{s} a_{ij} c_{j}^{k-1} = \frac{\bar{c}_i^k}{k}  \qquad i = 1 , \, ... , \, \sigma , \; k = 1 , \, ... , \, s ,
```
where $c$ are Gauß-Legendre nodes with $s$ stages and $\bar{c}$ are Gauß-Lobatto nodes with $\sigma$ stages.
"""
function get_lobatto_glrk_coefficients(s, σ=s+1, T=Float64)
    if σ == 1
        @error "Lobatto III coefficients for one stage are not defined."
    end

    function get_glrk_nodes(s)
        if s == 1
            c = [1//2]
        elseif s == 2
            c = [1/2-√3/6, 1/2+√3/6]
        elseif s == 3
            c = [1/2-√15/10,  1/2,        1/2+√15/10 ]
        elseif s == 4
            c = [
                    1/2-√(3/7+2*√30/35)/2,
                    1/2-√(3/7-2*√30/35)/2,
                    1/2+√(3/7-2*√30/35)/2,
                    1/2+√(3/7+2*√30/35)/2
                ]
        else
            c = nodes(GaussLegendreQuadrature(s))
        end

        return big.(c)
    end

    c = get_glrk_nodes(s)
    b̄ = get_lobatto_weights(σ)
    c̄ = get_lobatto_nodes(σ)
    M = [ c[j]^(k-1) for k in 1:s, j in 1:s ]
    
    row(i) = begin
        r = [ c̄[i]^k / k for k in 1:s ]
        M \ r
    end
    
    ā = vcat([row(i)' for i in 1:σ]...)

    CoefficientsIRK{T}(:LobattoIIIGLRK, s^2, s, σ, ā, b̄, c̄)
end


function get_lobatto_d_vector(s, T=Float64; normalize=false)
    q = get_lobatto_nodes(s)
    l = LagrangeBasis(q)
    v = [deriv_basis(l, j, i) for i in 1:s, j in 1:s]
    w = LinearAlgebra.nullspace(v')[:,begin]
    normalize ? T.(LinearAlgebra.normalize(w) .* sign(w[begin])) : T.(w)
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
