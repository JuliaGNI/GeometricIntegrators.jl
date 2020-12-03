
import LinearAlgebra
import Polynomials

import GeometricIntegrators.Quadratures: GaussLegendreQuadrature


@doc raw"""
The s-stage Lobatto nodes are defined as the roots of the following polynomial of degree $s$:
```math
\frac{d^{s-2}}{dx^{s-2}} \big( (x - x^2)^{s-1} \big) .
```
"""
function get_lobatto_nodes(s, T=BigFloat)
    if s == 1
        throw(ErrorException("Lobatto nodes for one stage are not defined."))
    end

    D(k) = Polynomials.derivative(Polynomial(T[0, 1, -1])^(k-1), k-2)
    c = sort(T.(Polynomials.roots(D(s))))
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
function get_lobatto_weights(s, T=BigFloat)
    if s == 1
        throw(ErrorException("Lobatto weights for one stage are not defined."))
    end

    P(k,x) = Polynomials.derivative(Polynomial(T[-1, 0, 1])^k, k)(x) / factorial(k) / 2^k
    c = get_lobatto_nodes(s,T)
    b = [ 1 / ( s*(s-1) * P(s-1, 2c[i] - 1)^2 ) for i in 1:s ]
end

@doc raw"""
The Lobatto IIIA coefficients are implicitly given by the so-called simplifying assumption $C(s)$:
```math
\sum \limits_{j=1}^{s} a_{ij} c_{j}^{k-1} = \frac{c_i^k}{k}  \qquad i = 1 , \, ... , \, s , \; k = 1 , \, ... , \, s .
```
"""
function get_lobatto_a_coefficients(s, T=BigFloat)
    if s == 1
        throw(ErrorException("Lobatto IIIA coefficients for one stage are not defined."))
    end
    solve_simplifying_assumption_c(get_lobatto_nodes(s,T))
end

@doc raw"""
The Lobatto IIIB coefficients are implicitly given by the so-called simplifying assumption $D(s)$:
```math
\sum \limits_{i=1}^{s} b_i c_{i}^{k-1} a_{ij} = \frac{b_j}{k} ( 1 - c_j^k)  \qquad j = 1 , \, ... , \, s , \; k = 1 , \, ... , \, s .
```
"""
function get_lobatto_b_coefficients(s, T=BigFloat)
    if s == 1
        throw(ErrorException("Lobatto IIIB coefficients for one stage are not defined."))
    end

    solve_simplifying_assumption_d(get_lobatto_weights(s,T), get_lobatto_nodes(s,T))
end

@doc raw"""
The Lobatto IIIC coefficients are determined by setting $a_{i,1} = b_1$ and
solving the so-called simplifying assumption $C(s-1)$, given by
```math
\sum \limits_{j=1}^{s} a_{ij} c_{j}^{k-1} = \frac{c_i^k}{k}  \qquad i = 1 , \, ... , \, s , \; k = 1 , \, ... , \, s-1 ,
```
for $a_{i,j}$ with $i = 1, ..., s$ and $j = 2, ..., s$.
"""
function get_lobatto_c_coefficients(s, T=BigFloat)
    if s == 1
        throw(ErrorException("Lobatto IIIC coefficients for one stage are not defined."))
    end

    b = get_lobatto_weights(s,T)
    c = get_lobatto_nodes(s,T)
    M = [ c[j]^(k-1) for k in 1:s-1, j in 2:s ]
    
    row(i) = begin
        r = [ c[i]^k / T(k) - c[1]^(k-1) * b[1] for k in 1:s-1 ]
        M \ r
    end
    
    hcat(b[1] * ones(T,s), vcat([transpose(row(i)) for i in 1:s]...))
end

@doc raw"""
The Lobatto IIIC̄ coefficients are determined by setting $a_{i,s} = 0$ and
solving the so-called simplifying assumption $C(s-1)$, given by
```math
\sum \limits_{j=1}^{s} a_{ij} c_{j}^{k-1} = \frac{c_i^k}{k}  \qquad i = 1 , \, ... , \, s , \; k = 1 , \, ... , \, s-1 ,
```
for $a_{i,j}$ with $i = 1, ..., s$ and $j = 1, ..., s-1$.
"""
function get_lobatto_c̄_coefficients(s, T=BigFloat)
    if s == 1
        throw(ErrorException("Lobatto IIIC̄ coefficients for one stage are not defined."))
    end

    c = get_lobatto_nodes(s,T)
    M = [ c[j]^(k-1) for k in 1:s-1, j in 1:s-1 ]
    
    row(i) = begin
        r = [ c[i]^k / T(k) for k in 1:s-1 ]
        M \ r
    end
    
    hcat(vcat([transpose(row(i)) for i in 1:s]...), zeros(T,s))
end

get_lobatto_d_coefficients(s, T=BigFloat) = (get_lobatto_c_coefficients(s,T) .+ get_lobatto_c̄_coefficients(s,T)) ./ 2

get_lobatto_e_coefficients(s, T=BigFloat) = (get_lobatto_a_coefficients(s,T) .+ get_lobatto_b_coefficients(s,T)) ./ 2

function get_lobatto_f_coefficients(s, T=BigFloat)
    if s == 1
        throw(ErrorException("Lobatto IIIF coefficients for one stage are not defined."))
    end

    c = get_lobatto_nodes(s,T)
    M = [ 1 / T(k + j - 1) for k in 1:s, j in 1:s ]
    r = [ 1 / T(s) / T(s + k) for k in 1:s ]
    α = M \ r
    
    Vₛ = [ c[i]^(j-1) for i in 1:s, j in 1:s ]
    Aₛ = zeros(T, s, s)
    for i in 2:s
       Aₛ[i,i-1] = 1 / T(i-1)
    end
    Aₛ[:,s] = α
    
    Vₛ * Aₛ * inv(Vₛ)
end

function get_lobatto_g_coefficients(s, T=BigFloat)
    a = get_lobatto_f_coefficients(s,T)
    b = get_lobatto_weights(s,T)
    ā = get_symplectic_conjugate_coefficients(a,b)
    return (a .+ ā) ./ 2
end


"Lobatto IIIA coefficients with s stages"
function CoefficientsLobattoIIIA(s, T=Float64)
    CoefficientsRK(T, Symbol("LobattoIIIA($s)"), 2s-2, get_lobatto_a_coefficients(s), get_lobatto_weights(s), get_lobatto_nodes(s))
end

"Lobatto IIIB coefficients with s stages"
function CoefficientsLobattoIIIB(s, T=Float64)
    CoefficientsRK(T, Symbol("LobattoIIIB($s)"), 2s-2, get_lobatto_b_coefficients(s), get_lobatto_weights(s), get_lobatto_nodes(s))
end

"Lobatto IIIC coefficients with s stages"
function CoefficientsLobattoIIIC(s, T=Float64)
    CoefficientsRK(T, Symbol("LobattoIIIC($s)"), 2s-2, get_lobatto_c_coefficients(s), get_lobatto_weights(s), get_lobatto_nodes(s))
end

"Lobatto IIIC̄ coefficients with s stages"
function CoefficientsLobattoIIIC̄(s, T=Float64)
    CoefficientsRK(T, Symbol("LobattoIIIC̄($s)"), 2s-2, get_lobatto_c̄_coefficients(s), get_lobatto_weights(s), get_lobatto_nodes(s))
end

"Lobatto IIID coefficients with s stages"
function CoefficientsLobattoIIID(s, T=Float64)
    CoefficientsRK(T, Symbol("LobattoIIID($s)"), 2s-2, get_lobatto_d_coefficients(s), get_lobatto_weights(s), get_lobatto_nodes(s))
end

"Lobatto IIIE coefficients with s stages"
function CoefficientsLobattoIIIE(s, T=Float64)
    CoefficientsRK(T, Symbol("LobattoIIIE($s)"), 2s-2, get_lobatto_e_coefficients(s), get_lobatto_weights(s), get_lobatto_nodes(s))
end

"Lobatto IIIF coefficients with s stages"
function CoefficientsLobattoIIIF(s, T=Float64)
    CoefficientsRK(T, Symbol("LobattoIIIF($s)"), 2s, get_lobatto_f_coefficients(s), get_lobatto_weights(s), get_lobatto_nodes(s))
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
