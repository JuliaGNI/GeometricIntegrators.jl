
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

"Construct splitting integrator."
function IntegratorSplitting(equation::SODE{DT,TT}, tableau::ST, Δt::TT) where {DT, TT, ST <: AbstractTableauSplitting{TT}}
    @assert has_exact_solution(equation)
    IntegratorSplitting{DT, ndims(equation)}(get_solution_tuple(equation), get_splitting_coefficients(length(equation.q), tableau)..., Δt)
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
