
import Polynomials
import Polynomials: Polynomial
import Polynomials: Polynomial
import SpecialPolynomials: ShiftedLegendre


function _ShiftedLegendre(s, T=BigFloat)
    p = vcat([zero(T) for _ in 1:s], one(T))
    ShiftedLegendre(p)
end


@doc raw"""
The Gauss nodes are given by the roots of the shifted Legendre polynomial
$P_s (2x-1)$ with $s$ the number of stages.
"""
function get_glrk_nodes(s, T=BigFloat)
    sort(T.(Polynomials.roots(convert(Polynomial, _ShiftedLegendre(s,T)))))
end

@doc raw"""
The Gauss weights are given by the following integrals
```math
b_i = \bigg( \frac{dP}{dx} (c_i) \bigg)^{-2} \int \limits_0^1 \bigg( \frac{P(x)}{x - c_i} \bigg)^2 dx ,
```
where $P(x)$ denotes the shifted Legendre polynomial
$P(x) = P_s (2x-1)$ with $s$ the number of stages.
"""
function get_glrk_weights(s, T=BigFloat)
    c = get_glrk_nodes(s,T)
    P = convert(Polynomial, _ShiftedLegendre(s,T))
    D = Polynomials.derivative(P)
    
    inti(i) = begin
        I = Polynomials.integrate( ( P ÷ Polynomial(T[-c[i], 1]) )^2 )
        I(1) - I(0)
    end
    
    b = [ inti(i) / D(c[i])^2  for i in 1:s ]
end

@doc raw"""
The Gauss coefficients are implicitly given by the so-called simplifying assumption $C(s)$:
```math
\sum \limits_{j=1}^{s} a_{ij} c_{j}^{k-1} = \frac{c_i^k}{k}  \qquad i = 1 , \, ... , \, s , \; k = 1 , \, ... , \, s .
```
"""
function get_glrk_coefficients(s, T=BigFloat)
    solve_simplifying_assumption_c(get_glrk_nodes(s,T))
end


function CoefficientsGLRK(s::Int; T=Float64)
    CoefficientsRK(T, Symbol("GLRK($s)"), 2s, get_glrk_coefficients(s), get_glrk_weights(s), get_glrk_nodes(s))
end


function get_GLRK_ω_matrix(s)
    ω = zeros(s, s+1)
    g = CoefficientsGLRK(s)
    ω[s,s+1] = 1
    for i in 1:s-1
        for j in 1:s
            ω[i,j] = g.b[j] * g.c[j]^(i-1)
        end
    end
    return ω
end
