
abstract type SplittingCoefficients end

coefficients(problem::SODEProblem, splitting::SplittingCoefficients) = coefficients(length(solutions(problem).q), splitting)


function _splitting_coefficients(r::Int, a::AbstractVector{T}, b::AbstractVector{T}) where {T}
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
General splitting method for vector fields with two components A and B.

Integrator:
```math
\varphi_{\tau} = \varphi_{b_s \tau}^{B} \circ \varphi_{a_s \tau}^{A} \circ \dotsc \circ \varphi_{b_1 \tau}^{B} \circ \varphi_{a_1 \tau}^{A}
```
"""
struct SplittingCoefficientsGeneral{T} <: SplittingCoefficients
    @HeaderTableau

    a::Vector{T}
    b::Vector{T}

    function SplittingCoefficientsGeneral{T}(name, o, s, a, b) where {T}
        @assert s == length(a) == length(b)
        new(name, o, s, a, b)
    end
end

function SplittingCoefficientsGeneral(name, o, a::Vector{T}, b::Vector{T}) where {T}
    SplittingCoefficientsGeneral{T}(name, o, length(a), a, b)
end

function coefficients(nequs::Int, method::SplittingCoefficientsGeneral{T}) where {T}
    @assert nequs == 2
    @assert length(method.a) == length(method.b)

    s = length(method.a)
    f = zeros(Int, 2s)
    c = zeros(T,   2s)

    for i in eachindex(method.a, method.b)
        f[2*(i-1)+1] = 1
        c[2*(i-1)+1] = method.a[i]
        f[2*(i-1)+2] = 2
        c[2*(i-1)+2] = method.b[i]
    end

    # select all entries with non-vanishing coefficients
    y = c .!= 0

    return f[y], c[y]
end


@doc raw"""
Non-symmetric splitting method.
    See McLachlan, Quispel, 2003, Equ. (4.10).
    The methods A and B are the composition of all vector fields in the SODE
    and its adjoint, respectively.

Basic method: Lie composition
```math
\begin{aligned}
\varphi_{\tau}^{A} &= \varphi_{\tau}^{v_1} \circ \varphi_{\tau}^{v_2} \circ \dotsc \circ \varphi_{\tau}^{v_{r-1}} \circ \varphi_{\tau}^{v_r} \\
\varphi_{\tau}^{B} &= \varphi_{\tau}^{v_r} \circ \varphi_{\tau}^{v_{r-1}} \circ \dotsc \circ \varphi_{\tau}^{v_2} \circ \varphi_{\tau}^{v_1}
\end{aligned}
```

Integrator:
```math
\varphi_{\tau}^{NS} = \varphi_{b_s \tau}^{B} \circ \varphi_{a_s \tau}^{A} \circ \dotsc \circ \varphi_{b_1 \tau}^{B} \circ \varphi_{a_1 \tau}^{A}
```
"""
struct SplittingCoefficientsNonSymmetric{T} <: SplittingCoefficients
    @HeaderTableau

    a::Vector{T}
    b::Vector{T}

    function SplittingCoefficientsNonSymmetric{T}(name, o, s, a, b) where {T}
        @assert s == length(a) == length(b)
        new(name, o, s, a, b)
    end
end

function SplittingCoefficientsNonSymmetric(name, o, a::Vector{T}, b::Vector{T}) where {T}
    SplittingCoefficientsNonSymmetric{T}(name, o, length(a), a, b)
end

function coefficients(nequs::Int, method::SplittingCoefficientsNonSymmetric)
    # R = length(equation.v)
    # S = method.s
    #
    # f = zeros(Int, 2R*S)
    # c = zeros(TT,  2R*S)
    #
    # for i in 1:S
    #     for j in 1:R
    #         f[(2i-2)*R+j] = j
    #         c[(2i-2)*R+j] = method.a[i]
    #     end
    #     for j in 1:R
    #         f[(2i-1)*R+j] = R-j+1
    #         c[(2i-1)*R+j] = method.b[i]
    #     end
    # end

    _splitting_coefficients(nequs::Int, method.a, method.b)
end


@doc raw"""
Symmetric splitting method with general stages.
    See McLachlan, Quispel, 2003, Equ. (4.11).

Basic method: Lie composition
```math
\begin{aligned}
\varphi_{\tau}^{A} &= \varphi_{\tau}^{v_1} \circ \varphi_{\tau}^{v_2} \circ \dotsc \circ \varphi_{\tau}^{v_{r-1}} \circ \varphi_{\tau}^{v_r} \\
\varphi_{\tau}^{B} &= \varphi_{\tau}^{v_r} \circ \varphi_{\tau}^{v_{r-1}} \circ \dotsc \circ \varphi_{\tau}^{v_2} \circ \varphi_{\tau}^{v_1}
\end{aligned}
```

Integrator:
```math
\varphi_{\tau}^{GS} = \varphi_{a_1 \tau}^{A} \circ \varphi_{b_1 \tau}^{B} \circ \dotsc \circ \varphi_{b_1 \tau}^{B} \circ \varphi_{a_1 \tau}^{A}
```
"""
struct SplittingCoefficientsGS{T} <: SplittingCoefficients
    @HeaderTableau

    a::Vector{T}
    b::Vector{T}

    function SplittingCoefficientsGS{T}(name, o, s, a, b) where {T}
        @assert s == length(a) == length(b)
        new(name, o, s, a, b)
    end
end

function SplittingCoefficientsGS(name, o, a::Vector{T}, b::Vector{T}) where {T}
    SplittingCoefficientsGS{T}(name, o, length(a), a, b)
end

function coefficients(nequs::Int, method::SplittingCoefficientsGS) 
    # R = nequs
    # S = method.s
    #
    # f = zeros(Int, 2R*S)
    # c = zeros(TT,  2R*S)
    #
    # for i in 1:S
    #     for j in 1:R
    #         f[(2i-2)*R+j] = j
    #         c[(2i-2)*R+j] = method.a[i]
    #     end
    #     for j in R:-1:1
    #         f[(2i-1)*R+j] = R-j+1
    #         c[(2i-1)*R+j] = method.b[i]
    #     end
    # end

    f, c = _splitting_coefficients(nequs::Int, method.a, method.b)
    vcat(f, f[end:-1:1]), vcat(c, c[end:-1:1])
end


@doc raw"""
Symmetric splitting method with symmetric stages.
    See McLachlan, Quispel, 2003, Equ. (4.6).

Basic method: symmetric Strang composition
```math
\varphi_{\tau}^{A} = \varphi_{\tau/2}^{v_1} \circ \varphi_{\tau/2}^{v_2} \circ \dotsc \circ \varphi_{\tau/2}^{v_{r-1}} \circ \varphi_{\tau/2}^{v_r} \circ \varphi_{\tau/2}^{v_r} \circ \varphi_{\tau/2}^{v_{r-1}} \circ \dotsc \circ \varphi_{\tau/2}^{v_2} \circ \varphi_{\tau/2}^{v_1}
```

Integrator:
```math
\varphi_{\tau}^{SS} = \varphi_{a_1 \tau}^{A} \circ \varphi_{a_2 \tau}^{A} \circ \dotsc \circ \varphi_{a_s \tau}^{A} \circ \dotsc \circ \varphi_{a_2 \tau}^{A} \circ \varphi_{a_1 \tau}^{A}
```
"""
struct SplittingCoefficientsSS{T} <: SplittingCoefficients
    @HeaderTableau

    a::Vector{T}

    function SplittingCoefficientsSS{T}(name, o, s, a) where {T}
        @assert s == length(a)
        new(name, o, s, a)
    end
end

function SplittingCoefficientsSS(name, o, a::Vector{T}) where {T}
    SplittingCoefficientsSS{T}(name, o, length(a), a)
end

function coefficients(nequs::Int, method::SplittingCoefficientsSS{T}) where {T}
    r = nequs
    a = vcat(method.a, method.a[end-1:-1:1]) ./ 2
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
