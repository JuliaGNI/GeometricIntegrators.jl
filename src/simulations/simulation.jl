
struct Simulation{ET <: GeometricEquation, IT <: Integrator, ST <: Solution}
    equation::ET
    integrator::IT
    solution::ST
    ncycle::Int
    run_id::String
    filename::String
end

@inline equation(sim::Simulation) = sim.equation
@inline integrator(sim::Simulation) = sim.integrator
@inline solution(sim::Simulation) = sim.solution
@inline cycles(sim::Simulation) = 1:sim.ncycle
@inline GeometricBase.eachsample(sim::Simulation) = eachsample(solution(sim))


function Simulation(equ::ET, int::IT, sol::ST, run_id::String, filename::String) where {ET,IT,ST}
    @assert mod(sol.ntime, sol.nwrite) == 0
    ncycle = div(sol.ntime, sol.nwrite)
    Simulation{ET,IT,ST}(equ, int, sol, ncycle, run_id, filename)
end

function Simulation(equ::GeometricEquation, int::Integrator, Δt, run_id, filename, ntime; nsave=DEFAULT_NSAVE, nwrite=DEFAULT_NWRITE)
    Simulation(equ, int, Solution(equ, Δt, ntime; nsave=nsave, nwrite=nwrite), run_id, filename)
end

function Simulation(equ::GeometricEquation, tableau::Union{AbstractTableau, Tableau}, Δt, run_id, filename, ntime; kwargs...)
    Simulation(equ, Integrator(equ, tableau, Δt), Δt, run_id, filename, ntime; kwargs...)
end

function Simulation(equ::GeometricEquation, integrator, tableau::Union{AbstractTableau, Tableau}, Δt, run_id, filename, ntime; kwargs...)
    Simulation(equ, integrator(equ, tableau, Δt), Δt, run_id, filename, ntime; kwargs...)
end


function run!(sim::Simulation)

    println("Running ", sim.run_id, "...")

    h5io = SolutionHDF5(sim.filename, solution(sim))

    # create atomic solution
    asol = AtomicSolution(solution(sim), integrator(sim))

    try
        # loop over integration cycles showing progress bar
        @showprogress 5 for c in cycles(sim)
            # loop over samples
            for m in eachsample(sim)
                # get cache from solution
                get_initial_conditions!(solution(sim), asol, m)

                # initilize integrator
                initialize!(integrator(sim), asol)

                # loop over time steps
                for n in eachtimestep(solution(sim))
                    integrate!(integrator(sim), solution(sim), asol, m, n)
                end
            end

            save(h5io, solution(sim))
            c == sim.ncycle || reset!(solution(sim))
        end
    catch ex
        if isa(ex, DomainError)
            @warn("DOMAIN ERROR")
        elseif isa(ex, ErrorException)
            @warn("Simulation exited early.")
            @warn(ex.msg)
        else
            close(h5io)
            throw(ex)
        end
    end

    close(h5io)

    return solution(sim)
end
