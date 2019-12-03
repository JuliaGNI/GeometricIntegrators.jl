
abstract type AbstractTableauSplitting{T <: Real} <: AbstractTableau{T} end


"Tableau for non-symmetric splitting methods."
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


"Tableau for symmetric splitting methods with general stages."
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


"Tableau for symmetric splitting methods with symmetric stages."
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


"Explicit Runge-Kutta integrator cache."
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
struct IntegratorSplitting{DT, TT, VT, ST <: AbstractTableauSplitting, FT, CT, N, D} <: DeterministicIntegrator{DT,TT}
    equation::SODE{DT,TT,VT,N}
    tableau::ST
    f::FT
    c::CT
    Δt::TT

    cache::IntegratorCacheSplitting{DT,TT,D}

    function IntegratorSplitting(equation::SODE{DT,TT,VT,N}, tableau::ST, f::Vector{Int}, c::Vector{TT}, Δt::TT) where {DT, TT, VT, ST <: AbstractTableauSplitting, N}
        @assert length(f) == length(c)
        ft = Tuple(f)
        ct = Tuple(c)
        D  = equation.d
        cache = IntegratorCacheSplitting{DT,TT,D}()
        new{DT,TT,VT,ST,typeof(ft),typeof(ct),N,D}(equation, tableau, ft, ct, Δt, cache)
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
function IntegratorSplitting(equation::SODE{DT,TT,VT}, tableau::ST, Δt::TT) where {DT, TT, VT, ST <: TableauSplittingNS{TT}}

    # basic method: Lie composition
    # \varphi_{\tau,A} = \varphi_{\tau,v_1} \circ \varphi_{\tau,v_2} \circ \hdots \varphi_{\tau,v_{r-1}} \circ \varphi_{\tau,v_r}
    # \varphi_{\tau,B} = \varphi_{\tau,v_r} \circ \varphi_{\tau,v_{r-1}} \circ \hdots \varphi_{\tau,v_2} \circ \varphi_{\tau,v_1}

    # integrator:
    # \varphi_{NS} = \varphi_{b_s \tau, B} \circ \varphi_{a_s \tau, A} \circ \hdots \circ \varphi_{b_1 \tau, B} \circ \varphi_{a_1 \tau, A}

    f, c = get_splitting_coefficients(length(equation.v), tableau.a, tableau.b)

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

    IntegratorSplitting(equation, tableau, f, c, Δt)
end


"Construct splitting integrator for symmetric splitting tableau with general stages."
function IntegratorSplitting(equation::SODE{DT,TT,VT}, tableau::ST, Δt::TT) where {DT, TT, VT, ST <: TableauSplittingGS{TT}}

    # basic method: Lie composition
    # \varphi_{\tau,A} = \varphi_{\tau,v_1} \circ \varphi_{\tau,v_2} \circ \hdots \varphi_{\tau,v_{r-1}} \circ \varphi_{\tau,v_r}
    # \varphi_{\tau,B} = \varphi_{\tau,v_r} \circ \varphi_{\tau,v_{r-1}} \circ \hdots \varphi_{\tau,v_2} \circ \varphi_{\tau,v_1}

    # integrator:
    # \varphi_{GS} = \varphi_{a_1 \tau, A} \circ \varphi_{b_1 \tau, B} \circ \hdots \circ \varphi_{b_1 \tau, B} \circ \varphi_{a_1 \tau, A}

    f, c = get_splitting_coefficients(length(equation.v), tableau.a, tableau.b)

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

    IntegratorSplitting(equation, tableau, vcat(f, f[end:-1:1]), vcat(c, c[end:-1:1]), Δt)
end


"Construct splitting integrator for symmetric splitting tableau with symmetric stages."
function IntegratorSplitting(equation::SODE{DT,TT,VT}, tableau::ST, Δt::TT) where {DT, TT, VT, ST <: TableauSplittingSS{TT}}

    # basic method: symmetric Strang composition
    # \varphi_{\tau,A} = \varphi_{\tau/2,v_1} \circ \varphi_{\tau/2,v_2} \circ \hdots \varphi_{\tau/2,v_{r-1}} \circ \varphi_{\tau/2,v_r}
    #              \circ \varphi_{\tau/2,v_r} \circ \varphi_{\tau/2,v_{r-1}} \circ \hdots \varphi_{\tau/2,v_2} \circ \varphi_{\tau/2,v_1}

    # integrator:
    # \varphi_{SS} = \varphi_{a_1 \tau, A} \circ \varphi_{a_2 \tau, A} \hdots \circ \varphi_{a_s \tau, A} \circ \hdots \circ \varphi_{a_2 \tau, A} \circ \varphi_{a_1 \tau, A}

    r = length(equation.v)
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

    IntegratorSplitting(equation, tableau, f, c, Δt)
end


equation(int::IntegratorSplitting) = int.equation
timestep(int::IntegratorSplitting) = int.Δt


"Integrate ODE with splitting integrator."
function integrate_step!(int::IntegratorSplitting{DT,TT,FT}, sol::AtomisticSolutionODE{DT,TT}) where {DT,TT,FT}
    local tᵢ::TT

    # reset cache
    reset!(sol, timestep(int))

    # compute internal stages
    for i in eachindex(int.f, int.c)
        if int.c[i] ≠ zero(TT)
            tᵢ = sol.t̅ + timestep(int) * int.c[i]
            int.equation.v[int.f[i]](tᵢ, sol.q, int.cache.q̃, int.c[i] * timestep(int))
            sol.q .= int.cache.q̃
        end
    end
end
