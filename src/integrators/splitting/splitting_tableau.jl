
abstract type AbstractTableauSplitting{T <: Real} <: AbstractTableau{T} end

function get_splitting_coefficients(r, a::Vector{T}, b::Vector{T}) where {T}
    @assert length(a) == length(b)

    s = length(a)
    f = zeros(Int, 2r*s)
    c = zeros(T,   2r*s)

    for i in 1:s
        for j in 1:r
            f[(2i-2)*r+j] = j
            c[(2i-2)*r+j] = a[i]
        end
        for j in 1:r
            f[(2i-1)*r+j] = r-j+1
            c[(2i-1)*r+j] = b[i]
        end
    end

    return f, c
end


@doc raw"""
Tableau for general splitting methods for vector fields with two components A and B.

Integrator:
```math
\varphi = \varphi_{b_s \tau, B} \circ \varphi_{a_s \tau, A} \circ \dotsc \circ \varphi_{b_1 \tau, B} \circ \varphi_{a_1 \tau, A}
```
"""
struct TableauSplitting{T} <: AbstractTableauSplitting{T}
    @HeaderTableau

    a::Vector{T}
    b::Vector{T}

    function TableauSplitting{T}(name, o, s, a, b) where {T}
        @assert s == length(a) == length(b)
        new(name, o, s, a, b)
    end
end

function TableauSplitting(name, o, a::Vector{T}, b::Vector{T}) where {T}
    TableauSplitting{T}(name, o, length(a), a, b)
end

function get_splitting_coefficients(nequs, tableau::TableauSplitting{T}) where {T}
    @assert nequs == 2
    @assert length(tableau.a) == length(tableau.b)

    s = length(tableau.a)
    f = zeros(Int, 2s)
    c = zeros(T,   2s)

    for i in eachindex(tableau.a, tableau.b)
        f[2*(i-1)+1] = 1
        c[2*(i-1)+1] = tableau.a[i]
        f[2*(i-1)+2] = 2
        c[2*(i-1)+2] = tableau.b[i]
    end

    # select all entries with non-vanishing coefficients
    y = c .!= 0

    return f[y], c[y]
end


@doc raw"""
Tableau for non-symmetric splitting methods.
    See McLachlan, Quispel, 2003, Equ. (4.10).
    The methods A and B are the composition of all vector fields in the SODE
    and its adjoint, respectively.

Basic method: Lie composition
```math
\begin{aligned}
\varphi_{\tau,A} &= \varphi_{\tau,v_1} \circ \varphi_{\tau,v_2} \circ \dotsc \circ \varphi_{\tau,v_{r-1}} \circ \varphi_{\tau,v_r} \\
\varphi_{\tau,B} &= \varphi_{\tau,v_r} \circ \varphi_{\tau,v_{r-1}} \circ \dotsc \circ \varphi_{\tau,v_2} \circ \varphi_{\tau,v_1}
\end{aligned}
```

Integrator:
```math
\varphi_{NS} = \varphi_{b_s \tau, B} \circ \varphi_{a_s \tau, A} \circ \dotsc \circ \varphi_{b_1 \tau, B} \circ \varphi_{a_1 \tau, A}
```
"""
struct TableauSplittingNS{T} <: AbstractTableauSplitting{T}
    @HeaderTableau

    a::Vector{T}
    b::Vector{T}

    function TableauSplittingNS{T}(name, o, s, a, b) where {T}
        @assert s == length(a) == length(b)
        new(name, o, s, a, b)
    end
end

function TableauSplittingNS(name, o, a::Vector{T}, b::Vector{T}) where {T}
    TableauSplittingNS{T}(name, o, length(a), a, b)
end

function get_splitting_coefficients(nequs, tableau::TableauSplittingNS)
    # R = length(equation.v)
    # S = tableau.s
    #
    # f = zeros(Int, 2R*S)
    # c = zeros(TT,  2R*S)
    #
    # for i in 1:S
    #     for j in 1:R
    #         f[(2i-2)*R+j] = j
    #         c[(2i-2)*R+j] = tableau.a[i]
    #     end
    #     for j in 1:R
    #         f[(2i-1)*R+j] = R-j+1
    #         c[(2i-1)*R+j] = tableau.b[i]
    #     end
    # end

    get_splitting_coefficients(nequs, tableau.a, tableau.b)
end


@doc raw"""
Tableau for symmetric splitting methods with general stages.
    See McLachlan, Quispel, 2003, Equ. (4.11).

Basic method: Lie composition
```math
\begin{aligned}
\varphi_{\tau,A} &= \varphi_{\tau,v_1} \circ \varphi_{\tau,v_2} \circ \dotsc \circ \varphi_{\tau,v_{r-1}} \circ \varphi_{\tau,v_r} \\
\varphi_{\tau,B} &= \varphi_{\tau,v_r} \circ \varphi_{\tau,v_{r-1}} \circ \dotsc \circ \varphi_{\tau,v_2} \circ \varphi_{\tau,v_1}
\end{aligned}
```

Integrator:
```math
\varphi_{GS} = \varphi_{a_1 \tau, A} \circ \varphi_{b_1 \tau, B} \circ \dotsc \circ \varphi_{b_1 \tau, B} \circ \varphi_{a_1 \tau, A}
```
"""
struct TableauSplittingGS{T} <: AbstractTableauSplitting{T}
    @HeaderTableau

    a::Vector{T}
    b::Vector{T}

    function TableauSplittingGS{T}(name, o, s, a, b) where {T}
        @assert s == length(a) == length(b)
        new(name, o, s, a, b)
    end
end

function TableauSplittingGS(name, o, a::Vector{T}, b::Vector{T}) where {T}
    TableauSplittingGS{T}(name, o, length(a), a, b)
end

function get_splitting_coefficients(nequs, tableau::TableauSplittingGS) 
    # R = nequs
    # S = tableau.s
    #
    # f = zeros(Int, 2R*S)
    # c = zeros(TT,  2R*S)
    #
    # for i in 1:S
    #     for j in 1:R
    #         f[(2i-2)*R+j] = j
    #         c[(2i-2)*R+j] = tableau.a[i]
    #     end
    #     for j in R:-1:1
    #         f[(2i-1)*R+j] = R-j+1
    #         c[(2i-1)*R+j] = tableau.b[i]
    #     end
    # end

    f, c = get_splitting_coefficients(nequs, tableau.a, tableau.b)
    vcat(f, f[end:-1:1]), vcat(c, c[end:-1:1])
end


@doc raw"""
Tableau for symmetric splitting methods with symmetric stages.
    See McLachlan, Quispel, 2003, Equ. (4.6).

Basic method: symmetric Strang composition
```math
\varphi_{\tau,A} = \varphi_{\tau/2,v_1} \circ \varphi_{\tau/2,v_2} \circ \dotsc \circ \varphi_{\tau/2,v_{r-1}} \circ \varphi_{\tau/2,v_r}
                \circ \varphi_{\tau/2,v_r} \circ \varphi_{\tau/2,v_{r-1}} \circ \dotsc \circ \varphi_{\tau/2,v_2} \circ \varphi_{\tau/2,v_1}
```

Integrator:
```math
\varphi_{SS} = \varphi_{a_1 \tau, A} \circ \varphi_{a_2 \tau, A} \dotsc \circ \varphi_{a_s \tau, A} \circ \dotsc \circ \varphi_{a_2 \tau, A} \circ \varphi_{a_1 \tau, A}
```
"""
struct TableauSplittingSS{T} <: AbstractTableauSplitting{T}
    @HeaderTableau

    a::Vector{T}

    function TableauSplittingSS{T}(name, o, s, a) where {T}
        @assert s == length(a)
        new(name, o, s, a)
    end
end

function TableauSplittingSS(name, o, a::Vector{T}) where {T}
    TableauSplittingSS{T}(name, o, length(a), a)
end

function get_splitting_coefficients(nequs, tableau::TableauSplittingSS{T}) where {T}
    r = nequs
    a = vcat(tableau.a, tableau.a[end-1:-1:1]) ./ 2
    s = length(a)

    f = zeros(Int, 2r*s)
    c = zeros(T,   2r*s)

    for i in 1:s
        for j in 1:r
            f[(2i-2)*r+j] = j
            c[(2i-2)*r+j] = a[i]
            f[(2i-1)*r+j] = r-j+1
            c[(2i-1)*r+j] = a[i]
        end
    end

    f, c
end
