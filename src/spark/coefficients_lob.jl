
@doc raw"""
The projective Lobatto-GLRK coefficients are implicitly given by
```math
\sum \limits_{j=1}^{s} a_{ij} c_{j}^{k-1} = \frac{\bar{c}_i^k}{k}  \qquad i = 1 , \, ... , \, \sigma , \; k = 1 , \, ... , \, s ,
```
where $c$ are Gauß-Legendre nodes with $s$ stages and $\bar{c}$ are Gauß-Lobatto nodes with $\sigma$ stages.
"""
function lobatto_gauss_coefficients(s, σ = s+1, T = Float64)
    if σ == 1
        @error "Lobatto III coefficients for one stage are not defined."
    end

    # Computed in BigFloat and converted to T only in the constructor below: the
    # Vandermonde system M is ill-conditioned, so the extra precision is what makes ā
    # accurate to Float64. Solving it in Float64 instead costs 1.1e-15 at s = 4 and
    # 3.4e-13 at s = 8. QuadratureRules defaults to Float64, unlike the RungeKutta
    # accessors these replaced, so the element type has to be passed explicitly. Its
    # nodes and weights are on [0,1] by default, which is the Runge-Kutta convention.
    c = gauss_legendre_nodes(BigFloat, s)
    b̄ = lobatto_legendre_weights(BigFloat, σ)
    c̄ = lobatto_legendre_nodes(BigFloat, σ)
    M = [c[j]^(k-1) for k in 1:s, j in 1:s]

    row(i) = begin
        r = [c̄[i]^k / k for k in 1:s]
        M \ r
    end

    ā = vcat([row(i)' for i in 1:σ]...)

    CoefficientsIRK{T}(:LobattoIIIGLRK, s^2, s, σ, ā, b̄, c̄)
end

function lobatto_ω_matrix(s)
    as = TableauLobattoIIIA(s).a[2:s, 1:s]
    es = zeros(s)
    es[s] = 1

    Q = vcat(hcat(as, zeros(s-1)), hcat(zeros(s)', 1))
    L = vcat(as, es')
    ω = inv(L) * Q

    return ω
end
