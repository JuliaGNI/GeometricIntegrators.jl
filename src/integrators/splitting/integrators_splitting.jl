
"Splitting integrator."
struct IntegratorSplitting{DT, TT, D, S, QT <: Tuple} <: DeterministicIntegrator{DT,TT}
    q::QT
    f::NTuple{S,Int64}
    c::NTuple{S,TT}
    Δt::TT

    function IntegratorSplitting{DT,D}(solutions::solType, f::Vector{Int}, c::Vector{TT}, Δt::TT) where {DT, TT, D, solType <: Tuple}
        @assert length(f) == length(c)
        ft = Tuple(f)
        ct = Tuple(c)
        new{DT,TT,D,length(f),solType}(solutions, ft, ct, Δt)
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

    # compute splitting steps
    for i in eachindex(int.f, int.c)
        if int.c[i] ≠ zero(TT)
            cᵢ = timestep(int) * int.c[i]
            tᵢ = sol.t̅ + cᵢ

            # reset atomic solution
            reset!(sol, cᵢ)

            # compute new solution
            int.q[int.f[i]](tᵢ, sol.q̅, sol.q, cᵢ)
        end
    end
end
