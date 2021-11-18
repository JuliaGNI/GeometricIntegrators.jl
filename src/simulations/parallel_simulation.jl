"""
   ParallelSimulation: collects all data structures that are necessary for a
   simulation run in parallel.

   WARNING: Many integrators are not thread-safe at this point!

"""
struct ParallelSimulation{ET <: Equation, IT <: Tuple, ST <: Solution}
    equation::ET
    integrators::IT
    solution::ST
    ncycle::Int
    run_id::String
    filename::String
end

@inline equation(sim::ParallelSimulation) = sim.equation
@inline integrators(sim::ParallelSimulation) = sim.integrators
@inline integrator(sim::ParallelSimulation, id) = sim.integrators[id]
@inline solution(sim::ParallelSimulation) = sim.solution
@inline cycles(sim::ParallelSimulation) = 1:sim.ncycle
@inline GeometricBase.eachsample(sim::ParallelSimulation) = eachsample(solution(sim))
@inline GeometricBase.eachsample(sim::ParallelSimulation, id::Int) = begin
    nthreads = Threads.nthreads()
    ns = nsamples(solution(sim))
    ns_thread = div(ns, nthreads)
    i1 = (id-1)*ns_thread + 1
    i2 =  id   *ns_thread
    id == nthreads ? (i1:ns) : (i1:i2)
end


function ParallelSimulation(equ::ET, ints::IT, sol::ST, run_id::String, filename::String) where {ET,IT,ST}
    @assert mod(sol.ntime, sol.nwrite) == 0
    ncycle = div(sol.ntime, sol.nwrite)
    ParallelSimulation{ET,IT,ST}(equ, ints, sol, ncycle, run_id, filename)
end

function ParallelSimulation(equ::Equation, ints::Tuple, Δt, run_id, filename, ntime; nsave=DEFAULT_NSAVE, nwrite=DEFAULT_NWRITE)
    ParallelSimulation(equ, ints, Solution(equ, Δt, ntime; nsave=nsave, nwrite=nwrite), run_id, filename)
end

function ParallelSimulation(equ::Equation, tableau::Union{AbstractTableau, Tableau}, Δt, run_id, filename, ntime; kwargs...)
    ints = Tuple(Integrator(equ, tableau, Δt) for i in 1:Threads.nthreads())
    ParallelSimulation(equ, ints, Δt, run_id, filename, ntime; kwargs...)
end

function ParallelSimulation(equ::Equation, integrator, tableau::Union{AbstractTableau, Tableau}, Δt, run_id, filename, ntime; kwargs...)
    ints = Tuple(integrator(equ, tableau, Δt) for i in 1:Threads.nthreads())
    ParallelSimulation(equ, ints, Δt, run_id, filename, ntime; kwargs...)
end


function run!(sim::ParallelSimulation)

    println("Running ", sim.run_id, "...")

    h5io = SolutionHDF5(sim.filename, solution(sim))

    try
        # loop over integration cycles showing progress bar
        @showprogress 5 for c in cycles(sim)
            # loop over samples
            Threads.@threads for m in eachsample(sim)
                id = Threads.threadid()

                # create atomic solution
                asol = AtomicSolution(solution(sim), integrator(sim, id))

                # get cache from solution
                get_initial_conditions!(solution(sim), asol, m)

                # initilize integrator
                initialize!(integrator(sim, id), asol)

                # loop over time steps
                for n in eachtimestep(solution(sim))
                    integrate!(integrator(sim, id), solution(sim), asol, m, n)
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
