
struct ParallelSimulation{ET <: Equation, IT <: Integrator, ST <: Solution}
    equation::ET
    integrator::IT
    solution::ST
    ncycle::Int
    run_id::String
    filename::String
end

equation(sim::ParallelSimulation) = sim.equation
integrator(sim::ParallelSimulation) = sim.integrator
solution(sim::ParallelSimulation) = sim.solution
cycles(sim::ParallelSimulation) = 1:sim.ncycle


function ParallelSimulation(equ::ET, int::IT, sol::ST, run_id::String, filename::String) where {ET,IT,ST}
    @assert mod(sol.ntime, sol.nwrite) == 0
    ncycle = div(sol.ntime, sol.nwrite)
    ParallelSimulation{ET,IT,ST}(equ, int, sol, ncycle, run_id, filename)
end

function ParallelSimulation(equ::Equation, int::Integrator, Δt, run_id, filename, ntime; nsave=1, nwrite=1)
    ParallelSimulation(equ, int, ParallelSolution(equ, Δt, ntime; nsave=nsave, nwrite=nwrite), run_id, filename)
end

function ParallelSimulation(equ::Equation, tableau::AbstractTableau, Δt, run_id, filename, ntime; kwargs...)
    ParallelSimulation(equ, Integrator(equ, tableau, Δt), Δt, run_id, filename, ntime; kwargs...)
end

function ParallelSimulation(equ::Equation, integrator, tableau::AbstractTableau, Δt, run_id, filename, ntime; kwargs...)
    ParallelSimulation(equ, integrator(equ, tableau, Δt), Δt, run_id, filename, ntime; kwargs...)
end


function run!(sim::ParallelSimulation)

    println("Running ", sim.run_id, "...")

    create_hdf5!(solution(sim), sim.filename)

    try
        # loop over integration cycles showing progress bar
        @showprogress 5 for c in cycles(sim)
            # loop over samples
            Threads.@threads for m in eachsample(solution(sim))
                # create atomic solution
                asol = AtomicSolution(equation(sim))

                # get cache from solution
                get_initial_conditions!(solution(sim), asol, m)

                # initilize integrator
                initialize!(integrator(sim), asol)

                # loop over time steps
                for n in eachtimestep(solution(sim))
                    integrate!(integrator(sim), solution(sim), asol, m, n)
                end
            end

            write_to_hdf5(solution(sim))
            reset!(solution(sim))
        end
    catch ex
        if isa(ex, DomainError)
            @warn("DOMAIN ERROR")
        elseif isa(ex, ErrorException)
            @warn("Simulation exited early.")
            @warn(ex.msg)
        else
            close(solution(sim))
            throw(ex)
        end
    end

    close(solution(sim))

    return solution(sim)
end
