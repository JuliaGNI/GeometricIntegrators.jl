
import LinearAlgebra

import RungeKutta: get_lobatto_nodes, get_lobatto_weights, get_gauss_nodes


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

    c = get_gauss_nodes(s)
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


function get_lobatto_ω_matrix(s)
    as = TableauLobattoIIIA(s).a[2:s,1:s]
    es = zeros(s)
    es[s] = 1

    Q = vcat( hcat(as, zeros(s-1)), hcat(zeros(s)', 1) )
    L = vcat( as, es' )
    ω = inv(L) * Q
    
    return ω
end
