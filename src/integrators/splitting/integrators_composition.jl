
"Composition integrator."
struct IntegratorComposition{DT, TT, D, S, IT <: Tuple} <: DeterministicIntegrator{DT,TT}
    ints::IT
    Δt::TT

    function IntegratorComposition{DT,D}(integrators::IT, Δt::TT) where {DT, TT, D, IT <: Tuple}
        new{DT, TT, D, length(integrators), IT}(integrators, Δt)
    end
end

"Construct composition integrator."
function IntegratorComposition(equation::SODE{DT,TT}, integrators::Tuple, tableau::ST, Δt::TT) where {DT, TT, ST <: AbstractTableauSplitting{TT}}
    D = ndims(equation)
    R = length(equation.v)

    f, c = get_splitting_coefficients(R, tableau)

    equations = get_function_tuple(equation)
    solutions = get_solution_tuple(equation)
    subints = ()

    # construct composition integrators
    for i in eachindex(f,c)
        if c[i] ≠ zero(TT)
            cᵢ = Δt * c[i]

            if integrators[f[i]] <: IntegratorExactODE
                @assert has_exact_solution(equation, f[i])
                subint = integrators[f[i]]{DT,D}(solutions[f[i]], cᵢ)
            else
                subint = integrators[f[i]]{DT,D}(equations[f[i]], cᵢ)
            end

            subints  = (subints..., subint)
        end
    end

    IntegratorComposition{DT, ndims(equation)}(subints, Δt)
end

function IntegratorComposition(equation::SODE, tableau::AbstractTableauSplitting, Δt::Number)
    @assert has_exact_solution(equation)
    integrators = Tuple(IntegratorExactODE for q in equation.q)
    IntegratorComposition(equation, integrators, tableau, Δt)
end


timestep(int::IntegratorComposition) = int.Δt


"Integrate ODE with Composition integrator."
function integrate_step!(int::IntegratorComposition{DT,TT}, sol::AtomicSolutionODE{DT,TT}) where {DT,TT}
    # compute composition steps
    for subint in int.ints
        integrate_step!(subint, sol)
    end
end
