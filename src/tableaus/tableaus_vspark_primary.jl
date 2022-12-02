@doc raw"""

Consider a symplectic pair of tableaus $(a^{1}, b^{1}, c^{1})$ and $(a^{3}, b^{3}, c^{3})$, i.e., satsifying $b^{1}_{i} b^{3}_{j} = b^{1}_{i} a^{3}_{ij} + b^{3}_{j} a^{1}_{ji}$, with an arbitrary number of stages $s$.
Use the same tableaus for $\tilde{a}^{1}$ and $\tilde{a}^{3}$, so that $\tilde{s} = s$, as well as

```math
\begin{aligned}
\begin{array}{c|cc}
& \tfrac{1}{2} b^{1} \\
& \vdots \\
& \tfrac{1}{2} b^{1} \\
\hline
a^{2} & \\
\end{array}
&&
\begin{array}{c|cc}
& \tfrac{1}{2} b^{3} \\
& \vdots \\
& \tfrac{1}{2} b^{3} \\
\hline
a^{4} & \\
\end{array}
&&
\begin{array}{c|cc}
& \tfrac{1}{2} b^{1} \\
c^{1} & \vdots \\
& \tfrac{1}{2} b^{1} \\
\hline
\tilde{a}^{2} & \tfrac{1}{2} (1 + R(\infty)) \, b^{1} \\
\end{array}
&&
\begin{array}{c|cc}
& \tfrac{1}{2} b^{3} \\
c^{3} & \vdots \\
& \tfrac{1}{2} b^{3} \\
\hline
\tilde{a}^{4} & \tfrac{1}{2} (1 + R(\infty)) \, b^{3} \\
\end{array}
\end{aligned}
```

Set $\omega = [0, ..., 0, 1]$ and
```math
\delta_{ij} = \begin{cases}
+1 & j = i , \\
-1 & j = \tilde{s} , \\
 0 & \text{else} ,
\end{cases}
```
so that $\Lambda_{1} = \Lambda_{2} = ... = \Lambda_{\tilde{s}}$.

This methods is constructed to satisfy the constraint on the projective stages, $\phi(\tilde{Q}_{n,i}, \tilde{P}_{n,i}) = 0$ for $i = 1, \, ..., \, \tilde{s}$.
Note, however, that it violates the symplecticity conditions $b^{1}_{i} b^{4}_{j} = b^{1}_{i} a^{4}_{ij} + b^{4}_{j} \tilde{a}^{1}_{ji}$ and $b^{2}_{i} b^{3}_{j} = b^{2}_{i} \tilde{a}^{3}_{ij} + b^{3}_{j} a^{2}_{ji}$.

"""
function TableauVSPARKInternalProjection(name, q::Tableau{T}, p::Tableau{T}, d=[]; R∞=1) where {T}

    @assert q.s == p.s

    s = q.s
    o = min(q.o, p.o)

    α_q = zeros(T, s, s)
    α_p = zeros(T, s, s)

    for i in 1:s
        α_q[i,:] .= q.b ./ 2
        α_p[i,:] .= p.b ./ 2
    end

    β_q = q.b .* (1 + R∞) / 2
    β_p = p.b .* (1 + R∞) / 2
    
    c_λ = q.c
    d_λ = q.b
    ω_λ = zeros(T, 1, s+1)
    δ_λ = zeros(T, s-1, s)
    
    ω_λ[1,s+1] = 1
    
    for i in 1:s-1
        δ_λ[i,i] = +1
        δ_λ[i,s] = -1
    end


    if length(d) == 0
        return TableauVSPARKprimary(name, o,
                            q.a, p.a, α_q, α_p,
                            q.a, p.a, α_q, α_p,
                            q.b, p.b, β_q, β_p,
                            q.c, p.c, c_λ, d_λ,
                            ω_λ, δ_λ)
    else
        @assert length(d) == q.s == p.s

        return TableauVSPARKprimary(name, o,
                            q.a, p.a, α_q, α_p,
                            q.a, p.a, α_q, α_p,
                            q.b, p.b, β_q, β_p,
                            q.c, p.c, c_λ, d_λ,
                            ω_λ, δ_λ, d)
    end

end

"Tableau for Gauss-Legendre method with s stages and symplectic projection."
function TableauVSPARKGLRKpInternal(s)
    glrk = TableauGauss(s)
    TableauVSPARKInternalProjection(Symbol("GLRK($s)pInternal"), glrk, glrk; R∞=(-1)^s)
end



@doc raw"""

Consider a symplectic pair of tableaus $(a^{1}, b^{1}, c^{1})$ and $(a^{3}, b^{3}, c^{3})$, i.e., satsifying $b^{1}_{i} b^{3}_{j} = b^{1}_{i} a^{3}_{ij} + b^{3}_{j} a^{1}_{ji}$, with an arbitrary number of stages $s$, and set

```math
\begin{aligned}
\begin{array}{c|cc}
& \tfrac{1}{2} b^{1} \\
& \vdots \\
& \tfrac{1}{2} b^{1} \\
\hline
a^{2} & \\
\end{array}
&&
\begin{array}{c|cc}
& \tfrac{1}{2} b^{4} \\
& \vdots \\
& \tfrac{1}{2} b^{4} \\
\hline
a^{4} & \\
\end{array}
&&
\begin{array}{c|cc}
& \tfrac{1}{2} b^{1} \\
c^{1} & \vdots \\
& \tfrac{1}{2} b^{1} \\
\hline
\tilde{a}^{2} & \tfrac{1}{2} (1 + R(\infty)) \, b^{1} \\
\end{array}
&&
\begin{array}{c|cc}
& \tfrac{1}{2} b^{3} \\
c^{3} & \vdots \\
& \tfrac{1}{2} b^{3} \\
\hline
\tilde{a}^{4} & \tfrac{1}{2} (1 + R(\infty)) \, b^{3} \\
\end{array}
\end{aligned}
```

Note that by this definition $\tilde{s} = s$.
The coefficients $\tilde{a}^{1}$ and $\tilde{a}^{3}$ are determined by the (modified) symplecticity conditions, specifically $a^{4}_{ij} = b^{3}_{j} ( b^{1}_{i} - \tilde{a}^{1}_{ji}) / b^{1}_{i}$ and $a^{2}_{ij} = b^{1}_{j} ( b^{3}_{i} - \tilde{a}^{3}_{ji} ) / b^{3}_{i}$, where $b^{2}$ has been replaced with $b^{1}$ and $b^{4}$ with $b^{3}$, respectively.
Set $\omega = [0, ..., 0, 1]$ and
```math
\delta_{ij} = \begin{cases}
+1 & j = i , \\
-1 & j = \tilde{s} , \\
 0 & \text{else} ,
\end{cases}
```
so that $\Lambda_{1} = \Lambda_{2} = ... = \Lambda_{\tilde{s}}$.

Note that this method satisfies the symplecticity conditions $b^{1}_{i} b^{4}_{j} = b^{1}_{i} a^{4}_{ij} + b^{4}_{j} \tilde{a}^{1}_{ji}$ and $b^{2}_{i} b^{3}_{j} = b^{2}_{i} \tilde{a}^{3}_{ij} + b^{3}_{j} a^{2}_{ji}$ only if $R(\infty) = 1$ due to the definitions of $b^{2}$ and $b^{4}$.
Moreover, it does usually not satisfy the constraint on the projective stages, $\phi(\tilde{Q}_{n,i}, \tilde{P}_{n,i}) = 0$ for $i = 1, \, ..., \, \tilde{s}$, exactly, but only approximately with bounded error, thus implying a residual in the symplecticity equation even if $R(\infty) = 1$.

"""
function TableauVSPARKModifiedInternalProjection(name, q, p, d=[]; R∞=1)
    T = eltype(q.a)

    @assert q.s == p.s

    s = q.s
    o = min(q.o, p.o)

    α_q = zeros(T, s, s)
    α_p = zeros(T, s, s)

    for i in 1:s
        α_q[i,:] .= q.b ./ 2
        α_p[i,:] .= p.b ./ 2
    end
    
    β_q = q.b .* (1 + R∞) / 2
    β_p = p.b .* (1 + R∞) / 2

#    a_q̃, a_p̃ = get_ã_vspark_primary(α_q, β_q, q.b, α_p, β_p, p.b)
    a_q̃, a_p̃ = get_ã_vspark_primary(α_q, q.b, q.b, α_p, p.b, p.b)

    c_λ = q.c
    d_λ = q.b
    ω_λ = zeros(T, 1, s+1)
    δ_λ = zeros(T, s-1, s)
    
    ω_λ[1,s+1] = 1
    
    for i in 1:s-1
        δ_λ[i,i] = +1
        δ_λ[i,s] = -1
    end

    if length(d) == 0
        return TableauVSPARKprimary(name, o,
                            q.a, p.a, α_q, α_p,
                            a_q̃, a_p̃, α_q, α_p,
                            q.b, p.b, β_q, β_p,
                            q.c, p.c, c_λ, d_λ,
                            ω_λ, δ_λ)
    else
        @assert length(d) == q.s == p.s

        return TableauVSPARKprimary(name, o,
                            q.a, p.a, α_q, α_p,
                            a_q̃, a_p̃, α_q, α_p,
                            q.b, p.b, β_q, β_p,
                            q.c, p.c, c_λ, d_λ,
                            ω_λ, δ_λ, d)
    end

end

"Tableau for Gauss-Legendre method with s stages and symplectic projection."
function TableauVSPARKGLRKpModifiedInternal(s)
    glrk = TableauGauss(s)
    TableauVSPARKModifiedInternalProjection(Symbol("GLRK($s)pModifiedInternal"), glrk, glrk; R∞=(-1)^s)
end



@doc raw"""

Consider a symplectic pair of tableaus $(a^{1}, b^{1}, c^{1})$ and $(a^{3}, b^{3}, c^{3})$, i.e., satsifying $b^{1}_{i} b^{3}_{j} = b^{1}_{i} a^{3}_{ij} + b^{3}_{j} a^{1}_{ji}$, with an arbitrary number of stages $s$.
For the projection, choose the tableau with $\tilde{s} = 1$ and $\rho = 0$, such that $\tilde{Q}_{n,1} = \tfrac{1}{2} ( q_{n} + q_{n+1})$, $\tilde{P}_{n,1} = \tfrac{1}{2} ( p_{n} + p_{n+1})$, i.e.,

```math
\begin{aligned}
\begin{array}{c|cc}
& \tfrac{1}{2} \\
& \vdots \\
& \tfrac{1}{2} \\
\hline
a^{2} & \\
\end{array}
&&
\begin{array}{c|cc}
& \tfrac{1}{2} \\
& \vdots \\
& \tfrac{1}{2} \\
\hline
a^{4} & \\
\end{array}
&&
\begin{array}{c|cc}
\tfrac{1}{2} & \tfrac{1}{2} \\
\hline
\tilde{a}^{2} & \tfrac{1}{2} (1 + R(\infty))\\
\end{array}
&&
\begin{array}{c|cc}
\tfrac{1}{2} & \tfrac{1}{2} \\
\hline
\tilde{a}^{4} & \tfrac{1}{2} (1 + R(\infty))\\
\end{array}
\end{aligned}
```

The coefficients $\tilde{a}^{1}$ and $\tilde{a}^{3}$ are determined by the symplecticity conditions, specifically $a^{4}_{ij} = b^{4}_{j} ( b^{1}_{i} - \tilde{a}^{1}_{ji}) / b^{1}_{i}$ and $a^{2}_{ij} = b^{2}_{j} ( b^{3}_{i} - \tilde{a}^{3}_{ji} ) / b^{3}_{i}$, and $\omega = [0, 1]$.

"""
function TableauVSPARKMidpointProjection(name, q::Tableau{T}, p::Tableau{T}, d=[]; R∞=1) where {T}
    @assert q.s == p.s
    s = q.s

    g = TableauGauss(1)
    α̃ = g.a
    β = Array(g.b)
    γ = g.c

    α = 0.5 * ones(T, s, 1)

    q_ã = zeros(T, 1, s)
    for i in 1:s
        q_ã[1,i] = q.b[i] / β[1] * ( β[1] - α[i,1] )
    end

    p_ã = zeros(T, 1, s)
    for i in 1:s
        p_ã[1,i] = p.b[i] / β[1] * ( β[1] - α[i,1] )
    end

    β .= (1 + R∞) / 2
    dλ = ones(T, 1)
    ω  = reshape(T[0  1], (1,2))
    δ  = zeros(T, 0, 1)

    if length(d) == 0
        return TableauVSPARKprimary(name, min(q.o, p.o),
                            q.a, p.a, α, α,
                            q_ã, p_ã, α̃, α̃,
                            q.b, p.b, β, β,
                            q.c, p.c, γ, dλ,
                            ω, δ)
    else
        @assert length(d) == q.s == p.s

        return TableauVSPARKprimary(name, min(q.o, p.o),
                            q.a, p.a, α, α,
                            q_ã, p_ã, α̃, α̃,
                            q.b, p.b, β, β,
                            q.c, p.c, γ, dλ,
                            ω, δ, d)
    end
end

"Tableau for Gauss-Lobatto IIIA-IIIB method with s stages and midpoint projection."
function TableauVSPARKLobattoIIIAIIIBpMidpoint(s)
    TableauVSPARKMidpointProjection(Symbol("LobattoIIIAIIIB($s)pMidpoint"), TableauLobattoIIIA(s), TableauLobattoIIIB(s), get_lobatto_nullvector(s); R∞=-1^(s+1))
end

"Tableau for Gauss-Lobatto IIIB-IIIA method with s stages and midpoint projection."
function TableauVSPARKLobattoIIIBIIIApMidpoint(s)
    TableauVSPARKMidpointProjection(Symbol("LobattoIIIBIIIA($s)pMidpoint"), TableauLobattoIIIB(s), TableauLobattoIIIA(s), get_lobatto_nullvector(s); R∞=-1^(s+1))
end

"Tableau for Gauss-Legendre method with s stages and midpoint projection."
function TableauVSPARKGLRKpMidpoint(s)
    glrk = TableauGauss(s)
    TableauVSPARKMidpointProjection(Symbol("GLRK($s)pMidpoint"), glrk, glrk; R∞=(-1)^s)
end



@doc raw"""

Consider a symplectic pair of tableaus $(a^{1}, b^{1}, c^{1})$ and $(a^{3}, b^{3}, c^{3})$, i.e., satsifying $b^{1}_{i} b^{3}_{j} = b^{1}_{i} a^{3}_{ij} + b^{3}_{j} a^{1}_{ji}$, with an arbitrary number of stages $s$.
For the projection, choose the tableau with $\tilde{s} = 1$ and $\rho = 0$, such that $\tilde{Q}_{n,1} = \tfrac{1}{2} ( q_{n} + q_{n+1})$, $\tilde{P}_{n,1} = \tfrac{1}{2} ( p_{n} + p_{n+1})$, i.e.,

```math
\begin{aligned}
\begin{array}{c|cc}
\tfrac{1}{2} & \tfrac{1}{2} b^{1} \\
\hline
\tilde{a}^{1} & \\
\end{array}
&&
\begin{array}{c|cc}
\tfrac{1}{2} & \tfrac{1}{2} \\
\hline
\tilde{a}^{2} & \tfrac{1}{2} ( 1 + R (\infty) ) \\
\end{array}
&&
\begin{array}{c|cc}
\tfrac{1}{2} & \tfrac{1}{2} b^{3} \\
\hline
\tilde{a}^{3} & \\
\end{array}
&&
\begin{array}{c|cc}
\tfrac{1}{2} & \tfrac{1}{2} \\
\hline
\tilde{a}^{4} & \tfrac{1}{2} ( 1 + R (\infty) ) \\
\end{array}
\end{aligned}
```

The coefficients $a^{2}$ and $a^{4}$ are determined by the symplecticity conditions, specifically $a^{4}_{ij} = b^{4}_{j} ( b^{1}_{i} - \tilde{a}^{1}_{ji}) / b^{1}_{i}$ and $a^{2}_{ij} = b^{2}_{j} ( b^{3}_{i} - \tilde{a}^{3}_{ji} ) / b^{3}_{i}$, and $\omega = [0, 1]$.

"""
function TableauVSPARKModifiedMidpointProjection(name, q::Tableau{T}, p::Tableau{T}, d=[]; R∞=1) where {T}
    @assert q.s == p.s
    s = q.s

    g = TableauGauss(1)
    α̃ = g.a
    β = Array(g.b)
    γ = g.c

    q_ã = reshape(q.b ./ 2, 1, s)
    p_ã = reshape(p.b ./ 2, 1, s)
    
    α_q = zeros(T, s, 1)
    α_p = zeros(T, s, 1)

    for i in 1:s
        α_q[i,1] = β[1] / p.b[i] * ( p.b[i] - p_ã[1,i] )
        α_p[i,1] = β[1] / q.b[i] * ( q.b[i] - q_ã[1,i] )
    end
    
    β .= (1 + R∞) / 2
    dλ = ones(T, 1)
    ω  = reshape(T[0  1], (1,2))
    δ  = zeros(T, 0, 1)

    if length(d) == 0
        return TableauVSPARKprimary(name, min(q.o, p.o),
                            q.a, p.a, α_q, α_p,
                            q_ã, p_ã, α̃, α̃,
                            q.b, p.b, β, β,
                            q.c, p.c, γ, dλ,
                            ω, δ)
    else
        @assert length(d) == q.s == p.s

        return TableauVSPARKprimary(name, min(q.o, p.o),
                            q.a, p.a, α_q, α_p,
                            q_ã, p_ã, α̃, α̃,
                            q.b, p.b, β, β,
                            q.c, p.c, γ, dλ,
                            ω, δ, d)
    end
end

"Tableau for Gauss-Lobatto IIIA-IIIB method with s stages and midpoint projection."
function TableauVSPARKLobattoIIIAIIIBpModifiedMidpoint(s)
    TableauVSPARKModifiedMidpointProjection(Symbol("LobattoIIIAIIIB($s)pModifiedMidpoint"), TableauLobattoIIIA(s), TableauLobattoIIIB(s), get_lobatto_nullvector(s); R∞=-1^(s+1))
end

"Tableau for Gauss-Lobatto IIIB-IIIA method with s stages and midpoint projection."
function TableauVSPARKLobattoIIIBIIIApModifiedMidpoint(s)
    TableauVSPARKModifiedMidpointProjection(Symbol("LobattoIIIBIIIA($s)pModifiedMidpoint"), TableauLobattoIIIB(s), TableauLobattoIIIA(s), get_lobatto_nullvector(s); R∞=-1^(s+1))
end

"Tableau for Gauss-Legendre method with s stages and midpoint projection."
function TableauVSPARKGLRKpModifiedMidpoint(s)
    glrk = TableauGauss(s)
    TableauVSPARKModifiedMidpointProjection(Symbol("GLRK($s)pModifiedMidpoint"), glrk, glrk; R∞=(-1)^s)
end



@doc raw"""

Consider a symplectic pair of tableaus $(a^{1}, b^{1}, c^{1})$ and $(a^{3}, b^{3}, c^{3})$, i.e., satsifying $b^{1}_{i} b^{3}_{j} = b^{1}_{i} a^{3}_{ij} + b^{3}_{j} a^{1}_{ji}$, with an arbitrary number of stages $s$.
For the projection, choose the tableau with $\tilde{s} = 2$ and $\rho = 1$, such that $\tilde{Q}_{n,1} = q_{n}$, $\tilde{Q}_{n,2} = q_{n+1}$, $\tilde{P}_{n,1} = p_{n}$, $\tilde{P}_{n,2} = p_{n+1}$, i.e.,

```math
\begin{aligned}
\begin{array}{c|cc}
0 & 0 \\
1 & b^{1} \\
\hline
\tilde{a}^{1} & \\
\end{array}
&&
\begin{array}{c|cc}
& 0 & 0 \\
& \tfrac{1}{2} & \tfrac{1}{2} \\
\hline
\tilde{a}^{2} & \tfrac{1}{2} & \tfrac{1}{2} \\
\end{array}
&&
\begin{array}{c|cc}
0 & 0 \\
1 & b^{3} \\
\hline
\tilde{a}^{3} & \\
\end{array}
&&
\begin{array}{c|cc}
& 0 & 0 \\
& \tfrac{1}{2} & \tfrac{1}{2} \\
\hline
\tilde{a}^{4} & \tfrac{1}{2} & \tfrac{1}{2} \\
\end{array}
\end{aligned}
```

The coefficients $a^{2}$ and $a^{4}$ are determined by the symplecticity conditions, specifically $a^{4}_{ij} = b^{4}_{j} ( b^{1}_{i} - \tilde{a}^{1}_{ji}) / b^{1}_{i}$ and $a^{2}_{ij} = b^{2}_{j} ( b^{3}_{i} - \tilde{a}^{3}_{ji} ) / b^{3}_{i}$.
Further choose $\omega = [1, 1, 0]$ and $\delta = [-1, R_{\infty}]$, so that $\tilde{\Lambda}_{n,1} = R_{\infty} \tilde{\Lambda}_{n,2}$ and
```math
\tilde{P}_{n,1} - \vartheta (\tilde{Q}_{n,1}) + R_{\infty} ( \tilde{P}_{n,2} - \vartheta (\tilde{Q}_{n,2}) ) = 0 .
```
Due to the particular choice of projective stages, this is equivalent to
```math
p_{n} - \vartheta (q_{n}) + R_{\infty} ( p_{n+1} - \vartheta (q_{n+1}) ) = 0 ,
```
so that the constraint $\phi(q_{n+1}, p_{n+1}) = 0$ is satisfied if $\phi(q_{n}, p_{n}) = 0$.
Note that the choice of $\tilde{a}^{2}$ and $\tilde{a}^{4}$ violates the symplecticity condition $b^{2}_{i} b^{4}_{j} = b^{2}_{i} \tilde{a}^{4}_{ij} + b^{4}_{j} \tilde{a}^{2}_{ji}$.

"""
function TableauVSPARKSymmetricProjection(name, q::Tableau{T}, p::Tableau{T}, d=[]; R∞=1) where {T}

    @assert q.s == p.s

    o = min(q.o, p.o)

    a_q = q.a
    a_p = p.a

    α_q = zeros(T, q.s, 2)
    α_q[:,1] .= 0.5

    α_p = zeros(T, p.s, 2)
    α_p[:,1] .= 0.5

    a_q̃ = vcat(zero(q.b)', q.b')
    a_p̃ = vcat(zero(p.b)', p.b')

    α_q̃ = [[0.0  0.0]
           [0.5  0.5]]
    α_p̃ = [[0.0  0.0]
           [0.5  0.5]]

    b_q = q.b
    b_p = p.b
    β_q = [0.5, 0.5]
    β_p = [0.5, 0.5]

    c_q = q.c
    c_p = p.c
    c_λ = [ 0.0, 1.0]
    d_λ = [ 0.5, 0.5*R∞]

    ω_λ = reshape(T[1  1  0], (1,3))
    δ_λ = reshape(T[-1   R∞], (1,2))


    if length(d) == 0
        return TableauVSPARKprimary(name, o,
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ, δ_λ)
    else
        @assert length(d) == q.s == p.s

        return TableauVSPARKprimary(name, o,
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ, δ_λ, d)
    end
end

"Tableau for Gauss-Lobatto IIIA-IIIB method with s stages and symmetric projection."
function TableauVSPARKLobattoIIIAIIIBpSymmetric(s)
    TableauVSPARKSymmetricProjection(Symbol("LobattoIIIAIIIB($s)pSymmetric"), TableauLobattoIIIA(s), TableauLobattoIIIB(s), get_lobatto_nullvector(s); R∞=-1^(s+1))
end

"Tableau for Gauss-Lobatto IIIB-IIIA method with s stages and symmetric projection."
function TableauVSPARKLobattoIIIBIIIApSymmetric(s)
    TableauVSPARKSymmetricProjection(Symbol("LobattoIIIBIIIA($s)pSymmetric"), TableauLobattoIIIB(s), TableauLobattoIIIA(s), get_lobatto_nullvector(s); R∞=-1^(s+1))
end

"Tableau for Gauss-Legendre method with s stages and symplectic projection."
function TableauVSPARKGLRKpSymmetric(s)
    glrk = TableauGauss(s)
    TableauVSPARKSymmetricProjection(Symbol("GLRK($s)pSymmetric"), glrk, glrk; R∞=(-1)^s)
end



@doc raw"""

Consider a symplectic pair of tableaus $(a^{1}, b^{1}, c^{1})$ and $(a^{3}, b^{3}, c^{3})$, i.e., satsifying $b^{1}_{i} b^{3}_{j} = b^{1}_{i} a^{3}_{ij} + b^{3}_{j} a^{1}_{ji}$, with an arbitrary number of stages $s$.
For the projection, choose the Lobatto-IIIA and IIIB tableaus with $\tilde{s} = 2$ stages for $(\tilde{a}^{4}, b^{4})$ and $(\tilde{a}^{2}, b^{2})$, respectively, and choose $\tilde{a}^{1}$ and $\tilde{a}^{3}$ such that the projective stages correspond to the initial condition and the solution, i.e.,

```math
\begin{aligned}
\begin{array}{c|cc}
0 & 0 \\
1 & b^{1} \\
\hline
\tilde{a}^{1} & \\
\end{array}
&&
\begin{array}{c|cc}
0 & \tfrac{1}{2} & 0 \\
1 & \tfrac{1}{2} & 0 \\
\hline
\tilde{a}^{2} & \\
\end{array}
&&
\begin{array}{c|cc}
0 & 0 \\
1 & b^{3} \\
\hline
\tilde{a}^{3} & \\
\end{array}
&&
\begin{array}{c|cc}
0 & 0 & 0 \\
1 & \tfrac{1}{2} & \tfrac{1}{2} \\
\hline
\tilde{a}^{4} & \\
\end{array}
\end{aligned}
```

and compute $a^{2}$ and $a^{4}$ by the symplecticity conditions, that is $a^{2}_{ij} = b^{2}_{j} ( b^{3}_{i} - \tilde{a}^{3}_{ji} ) / b^{3}_{i}$ and $a^{4}_{ij} = b^{4}_{j} ( b^{1}_{i} - \tilde{a}^{1}_{ji}) / b^{1}_{i}$.
Finally choose $\omega = [0, 0, 1]$ and $\delta = [-1, R_{\infty}]$, implying that $\rho = 1$.
By construction, this method satisfies all symplecticity conditions, but the constraint on the projection stages, $\phi(\tilde{Q}_{n,i}, \tilde{P}_{n,i}) = 0$ for $i = 1, \, ..., \, \tilde{s}$, is not satisfied exactly, but only approximately, although with bounded error.

"""
function TableauVSPARKLobattoIIIAIIIBProjection(name, q::Tableau{T}, p::Tableau{T}, d=[]; R∞=1) where {T}

    @assert q.s == p.s

    o = min(q.o, p.o)

    loba = TableauLobattoIIIA(2)
    lobb = TableauLobattoIIIB(2)

    a_q = q.a
    a_p = p.a
    b_q = q.b
    b_p = p.b

    a_q̃ = vcat(zero(q.b)', q.b')
    a_p̃ = vcat(zero(p.b)', p.b')

    α_q̃ = loba.a
    α_p̃ = lobb.a

    β_q = loba.b
    β_p = lobb.b

    α_q, α_p = get_α_vspark_primary(a_q̃, b_q, β_q, a_p̃, b_p, β_p)

    c_q = q.c
    c_p = p.c
    c_λ = [ 0.0, 1.0]
    d_λ = [ 0.5, 0.5]

    ω_λ = reshape(T[0  0  1], (1,3))
    δ_λ = reshape(T[-1   R∞], (1,2))


    if length(d) == 0
        return TableauVSPARKprimary(name, o,
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ, δ_λ)
    else
        @assert length(d) == q.s == p.s

        return TableauVSPARKprimary(name, o,
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ, δ_λ, d)
    end

end

"Tableau for Gauss-Lobatto IIIA-IIIB method with s stages and Lobatto-IIIA-IIIB projection."
function TableauVSPARKLobattoIIIAIIIBpLobattoIIIAIIIB(s)
    TableauVSPARKLobattoIIIAIIIBProjection(Symbol("LobattoIIIAIIIB($s)pLobattoIIIAIIIB(2)"), TableauLobattoIIIA(s), TableauLobattoIIIB(s), get_lobatto_nullvector(s); R∞=-1^(s+1))
end

"Tableau for Gauss-Lobatto IIIB-IIIA method with s stages and Lobatto-IIIA-IIIB projection."
function TableauVSPARKLobattoIIIBIIIApLobattoIIIAIIIB(s)
    TableauVSPARKLobattoIIIAIIIBProjection(Symbol("LobattoIIIBIIIA($s)pLobattoIIIAIIIB(2)"), TableauLobattoIIIB(s), TableauLobattoIIIA(s), get_lobatto_nullvector(s); R∞=-1^(s+1))
end

"Tableau for Gauss-Legendre method with s stages and Lobatto-IIIA-IIIB projection."
function TableauVSPARKGLRKpLobattoIIIAIIIB(s)
    glrk = TableauGauss(s)
    TableauVSPARKLobattoIIIAIIIBProjection(Symbol("GLRK($s)pLobattoIIIAIIIB(2)"), glrk, glrk; R∞=(-1)^s)
end



@doc raw"""

This methods is the same as `TableauVSPARKLobattoIIIAIIIBProjection`, except for using Lobatto-IIIA and IIIB tableaus with $\tilde{s} = 2$ stages for $(\tilde{a}^{2}, b^{2})$, and $(\tilde{a}^{4}, b^{4})$ respectively, instead of the other way around.

"""
function TableauVSPARKLobattoIIIBIIIAProjection(name, q::Tableau{T}, p::Tableau{T}, d=[]; R∞=1) where {T}

    @assert q.s == p.s

    o = min(q.o, p.o)

    loba = TableauLobattoIIIA(2)
    lobb = TableauLobattoIIIB(2)

    a_q = q.a
    a_p = p.a
    b_q = q.b
    b_p = p.b

    a_q̃ = vcat(zero(q.b)', q.b')
    a_p̃ = vcat(zero(p.b)', p.b')

    α_q̃ = lobb.a
    α_p̃ = loba.a

    β_q = lobb.b
    β_p = loba.b

    α_q, α_p = get_α_vspark_primary(a_q̃, b_q, β_q, a_p̃, b_p, β_p)

    c_q = q.c
    c_p = p.c
    c_λ = [ 0.0, 1.0]
    d_λ = [ 0.5, 0.5]

    ω_λ = reshape(T[0  0  1], (1,3))
    δ_λ = reshape(T[-1   R∞], (1,2))


    if length(d) == 0
        return TableauVSPARKprimary(name, o,
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ, δ_λ)
    else
        @assert length(d) == q.s == p.s

        return TableauVSPARKprimary(name, o,
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ, δ_λ, d)
    end

end

"Tableau for Gauss-Lobatto IIIA-IIIB method with s stages and Lobatto-IIIB-IIIA projection."
function TableauVSPARKLobattoIIIAIIIBpLobattoIIIBIIIA(s)
    TableauVSPARKLobattoIIIBIIIAProjection(Symbol("LobattoIIIAIIIB($s)pLobattoIIIBIIIA(2)"), TableauLobattoIIIA(s), TableauLobattoIIIB(s), get_lobatto_nullvector(s); R∞=-1^(s+1))
end

"Tableau for Gauss-Lobatto IIIB-IIIA method with s stages and Lobatto-IIIB-IIIA projection."
function TableauVSPARKLobattoIIIBIIIApLobattoIIIBIIIA(s)
    TableauVSPARKLobattoIIIBIIIAProjection(Symbol("LobattoIIIBIIIA($s)pLobattoIIIBIIIA(2)"), TableauLobattoIIIB(s), TableauLobattoIIIA(s), get_lobatto_nullvector(s); R∞=-1^(s+1))
end

"Tableau for Gauss-Legendre method with s stages and Lobatto-IIIB-IIIA projection."
function TableauVSPARKGLRKpLobattoIIIBIIIA(s)
    glrk = TableauGauss(s)
    TableauVSPARKLobattoIIIBIIIAProjection(Symbol("GLRK($s)pLobattoIIIBIIIA(2)"), glrk, glrk; R∞=(-1)^s)
end



@doc raw"""

Consider a symplectic pair of tableaus $(a^{1}, b^{1}, c^{1})$ and $(a^{3}, b^{3}, c^{3})$, i.e., satsifying $b^{1}_{i} b^{3}_{j} = b^{1}_{i} a^{3}_{ij} + b^{3}_{j} a^{1}_{ji}$, with an arbitrary number of stages $s$.
For the projection, choose the Lobatto-IIIA and IIIB tableaus with $\tilde{s} = 2$ stages for $(\tilde{a}^{4}, b^{4})$ and $(\tilde{a}^{2}, b^{2})$, respectively.

The coefficients $\tilde{a}^{1}$ and $\tilde{a}^{3}$ are determined by the relations
```math
\begin{aligned}
\sum \limits_{j=1}^{s} \tilde{a}^{1}_{ij} (c_{j}^{1})^{k-1} &= \frac{(c_{i}^{2})^k}{k} , \qquad &
\sum \limits_{j=1}^{s} \tilde{a}^{3}_{ij} (c_{j}^{3})^{k-1} &= \frac{(c_{i}^{4})^k}{k} , \qquad &
i &= 1 , \, ... , \, \tilde{s} , \qquad &
k &= 1 , \, ... , \, s .
\end{aligned}
```

The coefficients $a^{2}$ and $a^{4}$ by the symplecticity conditions, that is $a^{2}_{ij} = b^{2}_{j} ( b^{3}_{i} - \tilde{a}^{3}_{ji} ) / b^{3}_{i}$ and $a^{4}_{ij} = b^{4}_{j} ( b^{1}_{i} - \tilde{a}^{1}_{ji}) / b^{1}_{i}$.
Finally choose $\omega = [0, 0, 1]$ and $\delta = [-1, R_{\infty}]$, implying that $\rho = 1$.
By construction, this method satisfies all symplecticity conditions, but the constraint on the projection stages, $\phi(\tilde{Q}_{n,i}, \tilde{P}_{n,i}) = 0$ for $i = 1, \, ..., \, \tilde{s}$, is not satisfied exactly, but only approximately, although with bounded error.

"""
function TableauVSPARKModifiedLobattoIIIAIIIBProjection(name, q::Tableau{T}, p::Tableau{T}, d=[]; R∞=1) where {T}

    @assert q.s == p.s

    s = q.s
    o = min(q.o, p.o)

    loba = TableauLobattoIIIA(2)
    lobb = TableauLobattoIIIB(2)

    a_q = q.a
    a_p = p.a
    b_q = q.b
    b_p = p.b

    β_q = loba.b
    β_p = lobb.b

    a_q̃ = get_lobatto_glrk_coefficients(s, 2, T).a
    a_p̃ = get_lobatto_glrk_coefficients(s, 2, T).a

    α_q̃ = loba.a
    α_p̃ = lobb.a

    α_q, α_p = get_α_vspark_primary(a_q̃, b_q, β_q, a_p̃, b_p, β_p)
    
    c_q = q.c
    c_p = p.c
    c_λ = loba.c
    d_λ = loba.b

    ω_λ = reshape(T[0  0  1], (1,3))
    δ_λ = reshape(T[-1   R∞], (1,2))


    if length(d) == 0
        return TableauVSPARKprimary(name, o,
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ, δ_λ)
    else
        @assert length(d) == q.s == p.s

        return TableauVSPARKprimary(name, o,
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ, δ_λ, d)
    end

end

"Tableau for Gauss-Lobatto IIIA-IIIB method with s stages and Lobatto-IIIA-IIIB projection."
function TableauVSPARKLobattoIIIAIIIBpModifiedLobattoIIIAIIIB(s)
    TableauVSPARKModifiedLobattoIIIAIIIBProjection(Symbol("LobattoIIIAIIIB($s)pModifiedLobattoIIIAIIIB(2)"), TableauLobattoIIIA(s), TableauLobattoIIIB(s), get_lobatto_nullvector(s); R∞=-1^(s+1))
end

"Tableau for Gauss-Lobatto IIIB-IIIA method with s stages and Lobatto-IIIA-IIIB projection."
function TableauVSPARKLobattoIIIBIIIApModifiedLobattoIIIAIIIB(s)
    TableauVSPARKModifiedLobattoIIIAIIIBProjection(Symbol("LobattoIIIBIIIA($s)pModifiedLobattoIIIAIIIB(2)"), TableauLobattoIIIB(s), TableauLobattoIIIA(s), get_lobatto_nullvector(s); R∞=-1^(s+1))
end

"Tableau for Gauss-Legendre method with s stages and symplectic projection."
function TableauVSPARKGLRKpModifiedLobattoIIIAIIIB(s)
    glrk = TableauGauss(s)
    TableauVSPARKModifiedLobattoIIIAIIIBProjection(Symbol("GLRK($s)pModifiedLobattoIIIAIIIB(2)"), glrk, glrk; R∞=(-1)^s)
end



@doc raw"""

This methods is the same as `TableauVSPARKModifiedLobattoIIIAIIIBProjection`, except for using Lobatto-IIIA and IIIB tableaus with $\tilde{s} = 2$ stages for $(\tilde{a}^{2}, b^{2})$, and $(\tilde{a}^{4}, b^{4})$ respectively, instead of the other way around.

"""
function TableauVSPARKModifiedLobattoIIIBIIIAProjection(name, q::Tableau{T}, p::Tableau{T}, d=[]; R∞=1) where {T}

    @assert q.s == p.s

    s = q.s
    o = min(q.o, p.o)

    loba = TableauLobattoIIIA(2)
    lobb = TableauLobattoIIIB(2)

    a_q = q.a
    a_p = p.a
    b_q = q.b
    b_p = p.b

    β_q = loba.b
    β_p = lobb.b

    a_q̃ = get_lobatto_glrk_coefficients(s, 2, T).a
    a_p̃ = get_lobatto_glrk_coefficients(s, 2, T).a

    α_q̃ = lobb.a
    α_p̃ = loba.a

    α_q, α_p = get_α_vspark_primary(a_q̃, b_q, β_q, a_p̃, b_p, β_p)
    
    c_q = q.c
    c_p = p.c
    c_λ = loba.c
    d_λ = loba.b

    ω_λ = reshape(T[0  0  1], (1,3))
    δ_λ = reshape(T[-1   R∞], (1,2))


    if length(d) == 0
        return TableauVSPARKprimary(name, o,
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ, δ_λ)
    else
        @assert length(d) == q.s == p.s

        return TableauVSPARKprimary(name, o,
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ, δ_λ, d)
    end

end

"Tableau for Gauss-Lobatto IIIA-IIIB method with s stages and Lobatto-IIIB-IIIA projection."
function TableauVSPARKLobattoIIIAIIIBpModifiedLobattoIIIBIIIA(s)
    TableauVSPARKModifiedLobattoIIIBIIIAProjection(Symbol("LobattoIIIAIIIB($s)pModifiedLobattoIIIBIIIA(2)"), TableauLobattoIIIA(s), TableauLobattoIIIB(s), get_lobatto_nullvector(s); R∞=-1^(s+1))
end

"Tableau for Gauss-Lobatto IIIB-IIIA method with s stages and Lobatto-IIIB-IIIA projection."
function TableauVSPARKLobattoIIIBIIIApModifiedLobattoIIIBIIIA(s)
    TableauVSPARKModifiedLobattoIIIBIIIAProjection(Symbol("LobattoIIIBIIIA($s)pModifiedLobattoIIIBIIIA(2)"), TableauLobattoIIIB(s), TableauLobattoIIIA(s), get_lobatto_nullvector(s); R∞=-1^(s+1))
end

"Tableau for Gauss-Legendre method with s stages and symplectic projection."
function TableauVSPARKGLRKpModifiedLobattoIIIBIIIA(s)
    glrk = TableauGauss(s)
    TableauVSPARKModifiedLobattoIIIBIIIAProjection(Symbol("GLRK($s)pModifiedLobattoIIIBIIIA(2)"), glrk, glrk; R∞=(-1)^s)
end




function TableauVSPARKSymplecticProjection(name, q::Tableau{T}, p::Tableau{T}, d=[]; R∞=+1) where {T}

    la = TableauLobattoIIIA(2)
    lb = TableauLobattoIIIB(2)

    @assert q.s == p.s
    @assert q.b == p.b
    @assert q.c == p.c

    @assert la.s == lb.s
    @assert la.b == lb.b

    s = q.s
    σ = la.s
    ρ = 0

    a_q = q.a
    a_p = p.a
    b_q = q.b
    b_p = p.b

    β_q = [0.5, R∞*0.5]
    β_p = [0.5, R∞*0.5]

    α_q = zeros(T, s, 2)
    α_q[:,1] .= 0.5

    α_p = zeros(T, s, 2)
    α_p[:,1] .= 0.5

    a_q̃ = zeros(T, σ, s)
    a_p̃ = zeros(T, σ, s)
    for i in 1:σ
        for j in 1:s
            a_q̃[i,j] = b_q[j] / β_p[i] * (β_p[i] - α_p[j,i])
            a_p̃[i,j] = b_p[j] / β_q[i] * (β_q[i] - α_q[j,i])
        end
    end

    α_q̃ = la.a
    α_p̃ = lb.a

    c_q = q.c
    c_p = p.c
    c_λ = la.c
    d_λ = [0.5, 0.5]


    ω_λ = [0.5 0.5 0.0
           0.0 0.0 1.0]

    δ_λ = zeros(T, ρ, σ)


    if length(d) == 0
        return TableauVSPARKprimary(name, min(q.o, p.o),
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ, δ_λ)
    else
        @assert length(d) == s

        return TableauVSPARKprimary(name, min(q.o, p.o),
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ, δ_λ, d)
    end

end

function TableauVSPARKGLRKpSymplectic(s)
    glrk = TableauGauss(s)
    TableauVSPARKSymplecticProjection(Symbol("GLRK($s)pSymplectic"), glrk, glrk; R∞=(-1)^s)
end


function TableauVSPARKLobABCCD(s=2, T=Float64)
    lob  = TableauLobattoIIIC̄(s)
    loba = TableauLobattoIIIA(s)
    lobb = TableauLobattoIIIB(s)
    lobc = TableauLobattoIIIC(s)
    lobd = TableauLobattoIIID(s)
    
    a_q = lobd.a
    a_p = lobd.a
    
    α_q = lobb.a
    α_p = lobb.a

    a_q̃ = loba.a
    a_p̃ = loba.a
    
    α_q̃ = lob.a
    α_p̃ = lobc.a

    
    b_q = lobd.b
    b_p = lobd.b

    β_q̃ = lobd.b
    β_p̃ = lobd.b
    
    c_q = lobd.c
    c_p = lobd.c

    γ = lobd.c

    d = ones(T, s)
    ω = get_lobatto_ω_matrix(s)
    δ = zeros(T, 0, s)

    return TableauVSPARKprimary(Symbol("LobABCCD($s)"), 2*s-2,
                        a_q, a_p, α_q, α_p,
                        a_q̃, a_p̃, α_q̃, α_p̃,
                        b_q, b_p, β_q̃, β_p̃,
                        c_q, c_p, γ, d,
                        ω, δ)
end

function TableauVSPARKLobABCCE(s=2, T=Float64)
    lob  = TableauLobattoIIIC̄(s)
    loba = TableauLobattoIIIA(s)
    lobb = TableauLobattoIIIB(s)
    lobc = TableauLobattoIIIC(s)
    lobd = TableauLobattoIIID(s)
    lobe = TableauLobattoIIIE(s)
    
    a_q = lobe.a
    a_p = lobe.a
    
    α_q = lobb.a
    α_p = lobb.a

    a_q̃ = loba.a
    a_p̃ = loba.a
    
    α_q̃ = lob.a
    α_p̃ = lobc.a

    
    b_q = lobd.b
    b_p = lobd.b

    β_q̃ = lobd.b
    β_p̃ = lobd.b
    
    c_q = lobd.c
    c_p = lobd.c

    γ = lobd.c

    d = ones(T, s)
    ω = get_lobatto_ω_matrix(s)
    δ = zeros(T, 0, s)

    return TableauVSPARKprimary(Symbol("LobABCCE($s)"), 2*s-2,
                        a_q, a_p, α_q, α_p,
                        a_q̃, a_p̃, α_q̃, α_p̃,
                        b_q, b_p, β_q̃, β_p̃,
                        c_q, c_p, γ, d,
                        ω, δ)
end

function TableauVSPARKLobABED(s=2, T=Float64)
    lob  = TableauLobattoIIIC̄(s)
    loba = TableauLobattoIIIA(s)
    lobb = TableauLobattoIIIB(s)
    lobc = TableauLobattoIIIC(s)
    lobd = TableauLobattoIIID(s)
    lobe = TableauLobattoIIIE(s)

    a_q = lobe.a
    a_p = lobe.a
    
    α_q = lobb.a
    α_p = lobb.a

    a_q̃ = loba.a
    a_p̃ = loba.a
    
    α_q̃ = lobd.a
    α_p̃ = lobd.a

    
    b_q = lob.b
    b_p = lob.b

    β_q̃ = lob.b
    β_p̃ = lob.b
    
    c_q = lob.c
    c_p = lob.c

    γ   = lob.c

    d = ones(T, s)
    ω = get_lobatto_ω_matrix(s)
    δ = zeros(T, 0, s)

    return TableauVSPARKprimary(Symbol("LobABED($s)"), 2*s-2,
                        a_q, a_p, α_q, α_p,
                        a_q̃, a_p̃, α_q̃, α_p̃,
                        b_q, b_p, β_q̃, β_p̃,
                        c_q, c_p, γ, d,
                        ω, δ)
end

function TableauVSPARKLobABDE(s=2, T=Float64)
    lob  = TableauLobattoIIIC̄(s)
    loba = TableauLobattoIIIA(s)
    lobb = TableauLobattoIIIB(s)
    lobc = TableauLobattoIIIC(s)
    lobd = TableauLobattoIIID(s)
    lobe = TableauLobattoIIIE(s)

    a_q = lobd.a
    a_p = lobd.a
    
    α_q = lobb.a
    α_p = lobb.a

    a_q̃ = loba.a
    a_p̃ = loba.a
    
    α_q̃ = lobe.a
    α_p̃ = lobe.a

    
    b_q = lob.b
    b_p = lob.b

    β_q̃ = lob.b
    β_p̃ = lob.b
    
    c_q = lob.c
    c_p = lob.c

    γ   = lob.c

    d = ones(T, s)
    ω = get_lobatto_ω_matrix(s)
    δ = zeros(T, 0, s)

    return TableauVSPARKprimary(Symbol("LobABDE($s)"), 2*s-2,
                        a_q, a_p, α_q, α_p,
                        a_q̃, a_p̃, α_q̃, α_p̃,
                        b_q, b_p, β_q̃, β_p̃,
                        c_q, c_p, γ, d,
                        ω, δ)
end

function TableauVSPARKLobABD(s=2, T=Float64)
    lob  = TableauLobattoIIIC̄(s)
    loba = TableauLobattoIIIA(s)
    lobb = TableauLobattoIIIB(s)
    lobc = TableauLobattoIIIC(s)
    lobd = TableauLobattoIIID(s)
    lobe = TableauLobattoIIIE(s)


    a_q = lobd.a
    a_p = lobd.a
    
    α_q = lobb.a
    α_p = lobb.a

    a_q̃ = loba.a
    a_p̃ = loba.a
    
    α_q̃ = loba.a
    α_p̃ = lobb.a

    
    b_q = lob.b
    b_p = lob.b

    β_q̃ = lob.b
    β_p̃ = lob.b
    
    c_q = lob.c
    c_p = lob.c

    γ   = lob.c

    d = ones(T, s)
    ω = get_lobatto_ω_matrix(s)
    δ = zeros(T, 0, s)

    return TableauVSPARKprimary(Symbol("LobABD($s)"), 2*s-2,
                        a_q, a_p, α_q, α_p,
                        a_q̃, a_p̃, α_q̃, α_p̃,
                        b_q, b_p, β_q̃, β_p̃,
                        c_q, c_p, γ, d,
                        ω, δ)
end

function TableauVSPARKLobABE(s=2, T=Float64)
    lob  = TableauLobattoIIIC̄(s)
    loba = TableauLobattoIIIA(s)
    lobb = TableauLobattoIIIB(s)
    lobc = TableauLobattoIIIC(s)
    lobd = TableauLobattoIIID(s)
    lobe = TableauLobattoIIIE(s)


    a_q = lobe.a
    a_p = lobe.a
    
    α_q = lobb.a
    α_p = lobb.a

    a_q̃ = loba.a
    a_p̃ = loba.a
    
    α_q̃ = loba.a
    α_p̃ = lobb.a

    
    b_q = lob.b
    b_p = lob.b

    β_q̃ = lob.b
    β_p̃ = lob.b
    
    c_q = lob.c
    c_p = lob.c

    γ   = lob.c

    d = ones(T, s)
    ω = get_lobatto_ω_matrix(s)
    δ = zeros(T, 0, s)

    return TableauVSPARKprimary(Symbol("LobABE($s)"), 2*s-2,
                        a_q, a_p, α_q, α_p,
                        a_q̃, a_p̃, α_q̃, α_p̃,
                        b_q, b_p, β_q̃, β_p̃,
                        c_q, c_p, γ, d,
                        ω, δ)
end

function TableauVSPARKLobDE(s=2, T=Float64)
    lob  = TableauLobattoIIIC̄(s)
    loba = TableauLobattoIIIA(s)
    lobb = TableauLobattoIIIB(s)
    lobc = TableauLobattoIIIC(s)
    lobd = TableauLobattoIIID(s)
    lobe = TableauLobattoIIIE(s)


    a_q = lobd.a
    a_p = lobd.a
    
    α_q = lobe.a
    α_p = lobe.a

    a_q̃ = lobe.a
    a_p̃ = lobe.a
    
    α_q̃ = lobd.a
    α_p̃ = lobd.a

    
    b_q = lob.b
    b_p = lob.b

    β_q̃ = lob.b
    β_p̃ = lob.b
    
    c_q = lob.c
    c_p = lob.c

    γ   = lob.c

    d = ones(T, s)
    ω = get_lobatto_ω_matrix(s)
    δ = zeros(T, 0, s)

    return TableauVSPARKprimary(Symbol("LobDE($s)"), 2*s-2,
                        a_q, a_p, α_q, α_p,
                        a_q̃, a_p̃, α_q̃, α_p̃,
                        b_q, b_p, β_q̃, β_p̃,
                        c_q, c_p, γ, d,
                        ω, δ)
end

function TableauVSPARKLobED(s=2, T=Float64)
    lob  = TableauLobattoIIIC̄(s)
    loba = TableauLobattoIIIA(s)
    lobb = TableauLobattoIIIB(s)
    lobc = TableauLobattoIIIC(s)
    lobd = TableauLobattoIIID(s)
    lobe = TableauLobattoIIIE(s)


    a_q = lobe.a
    a_p = lobe.a
    
    α_q = lobd.a
    α_p = lobd.a

    a_q̃ = lobd.a
    a_p̃ = lobd.a
    
    α_q̃ = lobe.a
    α_p̃ = lobe.a

    
    b_q = lob.b
    b_p = lob.b

    β_q̃ = lob.b
    β_p̃ = lob.b
    
    c_q = lob.c
    c_p = lob.c

    γ   = lob.c

    d = ones(T, s)
    ω = get_lobatto_ω_matrix(s)
    δ = zeros(T, 0, s)

    return TableauVSPARKprimary(Symbol("LobED($s)"), 2*s-2,
                        a_q, a_p, α_q, α_p,
                        a_q̃, a_p̃, α_q̃, α_p̃,
                        b_q, b_p, β_q̃, β_p̃,
                        c_q, c_p, γ, d,
                        ω, δ)
end
