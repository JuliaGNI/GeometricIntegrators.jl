
struct Simulation{ET <: Equation, IT <: Integrator, ST <: Solution}
    equation::ET
    integrator::IT
    solution::ST
    run_id::String
    filename::String
end

function Simulation(equ::ET, int::IT, sol::ST, run_id::String, filename::String) where {ET,IT,ST}
    Simulation{ET,IT,TT}(equ, int, sol, run_id, filename)
end

function Simulation(equ::Equation, int::Integrator, Δt, run_id, filename, ntime, nsave)
    Simulation(equ, int, Solution(equ, Δt, ntime, nsave), run_id, filename)
end

function Simulation(equ::Equation, tableau::AbstractTableau, Δt, run_id, filename, ntime, nsave)
    Simulation(equ, Integrator(equ, tableau, Δt), Δt, run_id, filename, ntime, nsave)
end

function Simulation(equ::Equation, integrator, tableau::AbstractTableau, Δt, run_id, filename, ntime, nsave)
    Simulation(equ, integrator(equ, tableau, Δt), Δt, run_id, filename, ntime, nsave)
end


function run!(sim::Simulation)

    println("Running ", sim.run_id, "...")

    h5 = create_hdf5(sim.solution, sim.filename)

    try
        integrate!(sim.integrator, sim.solution)
    catch ex
        if isa(ex, DomainError)
            warn("DOMAIN ERROR")
        elseif isa(ex, ErrorException)
            warn("Simulation exited early.")
            warn(ex.msg)
        else
            throw(ex)
        end
    end

    write_to_hdf5(sim.solution, h5)
    close(h5)

    return sim.solution
end
