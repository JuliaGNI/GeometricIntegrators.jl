
abstract type AbstractTableauSplitting{T <: Real} <: AbstractTableau{T} end


"""
Tableau for non-symmetric splitting methods.
    See McLachlan, Quispel, 2003, Equ. (4.10).
    The methods A and B are the composition of all vector fields in the SODE
    and its adjoint, respectively.
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


"""
Tableau for symmetric splitting methods with general stages.
    See McLachlan, Quispel, 2003, Equ. (4.11).
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


"""
Tableau for symmetric splitting methods with symmetric stages.
    See McLachlan, Quispel, 2003, Equ. (4.6).
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


"Splitting integrator cache."
mutable struct IntegratorCacheSplitting{DT,TT,D} <: ODEIntegratorCache{DT,D}
    v::Vector{DT}
    q̃::Vector{DT}
    s̃::Vector{DT}

    function IntegratorCacheSplitting{DT,TT,D}() where {DT,TT,D}
        v = zeros(DT, D)
        q̃ = zeros(DT, D)
        s̃ = zeros(DT, D)
        new(v, q̃, s̃)
    end
end


"Splitting integrator."
struct IntegratorSplitting{DT, TT, D, S, QT <: Tuple} <: DeterministicIntegrator{DT,TT}
    q::QT
    f::NTuple{S,Int64}
    c::NTuple{S,TT}
    Δt::TT

    cache::IntegratorCacheSplitting{DT,TT,D}

    function IntegratorSplitting{DT,D}(solutions::solType, f::Vector{Int}, c::Vector{TT}, Δt::TT) where {DT, TT, D, solType <: Tuple}
        @assert length(f) == length(c)
        ft = Tuple(f)
        ct = Tuple(c)
        cache = IntegratorCacheSplitting{DT,TT,D}()
        new{DT,TT,D,length(f),solType}(solutions, ft, ct, Δt, cache)
    end
end


function get_splitting_coefficients(r, a::Vector{TT}, b::Vector{TT}) where {TT}
    @assert length(a) == length(b)

    s = length(a)
    f = zeros(Int, 2r*s)
    c = zeros(TT,  2r*s)

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


"Construct splitting integrator for non-symmetric splitting tableau with general stages."
function IntegratorSplitting(equation::SODE{DT,TT}, tableau::ST, Δt::TT) where {DT, TT, ST <: TableauSplittingNS{TT}}
    @assert has_exact_solution(equation)

    # basic method: Lie composition
    # \varphi_{\tau,A} = \varphi_{\tau,v_1} \circ \varphi_{\tau,v_2} \circ \hdots \varphi_{\tau,v_{r-1}} \circ \varphi_{\tau,v_r}
    # \varphi_{\tau,B} = \varphi_{\tau,v_r} \circ \varphi_{\tau,v_{r-1}} \circ \hdots \varphi_{\tau,v_2} \circ \varphi_{\tau,v_1}

    # integrator:
    # \varphi_{NS} = \varphi_{b_s \tau, B} \circ \varphi_{a_s \tau, A} \circ \hdots \circ \varphi_{b_1 \tau, B} \circ \varphi_{a_1 \tau, A}

    f, c = get_splitting_coefficients(length(equation.q), tableau.a, tableau.b)

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

    IntegratorSplitting{DT, ndims(equation)}(get_solution_tuple(equation), f, c, Δt)
end


"Construct splitting integrator for symmetric splitting tableau with general stages."
function IntegratorSplitting(equation::SODE{DT,TT}, tableau::ST, Δt::TT) where {DT, TT, ST <: TableauSplittingGS{TT}}
    @assert has_exact_solution(equation)

    # basic method: Lie composition
    # \varphi_{\tau,A} = \varphi_{\tau,v_1} \circ \varphi_{\tau,v_2} \circ \hdots \varphi_{\tau,v_{r-1}} \circ \varphi_{\tau,v_r}
    # \varphi_{\tau,B} = \varphi_{\tau,v_r} \circ \varphi_{\tau,v_{r-1}} \circ \hdots \varphi_{\tau,v_2} \circ \varphi_{\tau,v_1}

    # integrator:
    # \varphi_{GS} = \varphi_{a_1 \tau, A} \circ \varphi_{b_1 \tau, B} \circ \hdots \circ \varphi_{b_1 \tau, B} \circ \varphi_{a_1 \tau, A}

    f, c = get_splitting_coefficients(length(equation.q), tableau.a, tableau.b)

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
    #     for j in R:-1:1
    #         f[(2i-1)*R+j] = R-j+1
    #         c[(2i-1)*R+j] = tableau.b[i]
    #     end
    # end

    IntegratorSplitting{DT, ndims(equation)}(get_solution_tuple(equation), vcat(f, f[end:-1:1]), vcat(c, c[end:-1:1]), Δt)
end


"Construct splitting integrator for symmetric splitting tableau with symmetric stages."
function IntegratorSplitting(equation::SODE{DT,TT}, tableau::ST, Δt::TT) where {DT, TT, ST <: TableauSplittingSS{TT}}
    @assert has_exact_solution(equation)

    # basic method: symmetric Strang composition
    # \varphi_{\tau,A} = \varphi_{\tau/2,v_1} \circ \varphi_{\tau/2,v_2} \circ \hdots \varphi_{\tau/2,v_{r-1}} \circ \varphi_{\tau/2,v_r}
    #              \circ \varphi_{\tau/2,v_r} \circ \varphi_{\tau/2,v_{r-1}} \circ \hdots \varphi_{\tau/2,v_2} \circ \varphi_{\tau/2,v_1}

    # integrator:
    # \varphi_{SS} = \varphi_{a_1 \tau, A} \circ \varphi_{a_2 \tau, A} \hdots \circ \varphi_{a_s \tau, A} \circ \hdots \circ \varphi_{a_2 \tau, A} \circ \varphi_{a_1 \tau, A}

    r = length(equation.q)
    a = vcat(tableau.a, tableau.a[end-1:-1:1]) ./ 2
    s = length(a)

    f = zeros(Int, 2r*s)
    c = zeros(TT,  2r*s)

    for i in 1:s
        for j in 1:r
            f[(2i-2)*r+j] = j
            c[(2i-2)*r+j] = a[i]
            f[(2i-1)*r+j] = r-j+1
            c[(2i-1)*r+j] = a[i]
        end
    end

    IntegratorSplitting{DT, ndims(equation)}(get_solution_tuple(equation), f, c, Δt)
end


timestep(int::IntegratorSplitting) = int.Δt


"Integrate ODE with splitting integrator."
function integrate_step!(int::IntegratorSplitting{DT,TT}, sol::AtomicSolutionODE{DT,TT}) where {DT,TT}
    local cᵢ::TT
    local tᵢ::TT

    # reset atomic solution
    reset!(sol, timestep(int))

    # compute splitting steps
    for i in eachindex(int.f, int.c)
        if int.c[i] ≠ zero(TT)
            cᵢ = timestep(int) * int.c[i]
            tᵢ = sol.t̅ + cᵢ
            int.q[int.f[i]](tᵢ, sol.q, int.cache.q̃, cᵢ)
            sol.q .= int.cache.q̃
        end
    end
end
