
struct Simulation{ET <: Equation, IT <: Integrator, ST <: Solution}
    equation::ET
    integrator::IT
    solution::ST
    ncycle::Int
    run_id::String
    filename::String
end

function Simulation(equ::ET, int::IT, sol::ST, run_id::String, filename::String; ncycle=1) where {ET,IT,ST}
    Simulation{ET,IT,TT}(equ, int, sol, ncycle, run_id, filename)
end

function Simulation(equ::Equation, int::Integrator, Δt, run_id, filename, ntime, nsave; ncycle=1)
    Simulation(equ, int, Solution(equ, Δt, ntime, nsave), run_id, filename; ncycle=ncycle)
end

function Simulation(equ::Equation, tableau::AbstractTableau, Δt, run_id, filename, ntime, nsave; ncycle=1)
    Simulation(equ, Integrator(equ, tableau, Δt), Δt, run_id, filename, ntime, nsave; ncycle=ncycle)
end

function Simulation(equ::Equation, integrator, tableau::AbstractTableau, Δt, run_id, filename, ntime, nsave; ncycle=1)
    Simulation(equ, integrator(equ, tableau, Δt), Δt, run_id, filename, ntime, nsave; ncycle=ncycle)
end


function run!(sim::Simulation)

    println("Running ", sim.run_id, "...")

    h5 = create_hdf5(sim.solution, sim.filename)

    try
        @showprogress 5 for c in 1:sim.ncycle
            integrate!(sim.integrator, sim.solution)
            write_to_hdf5(sim.solution, h5)
            reset!(sim.solution)
        end
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

    close(h5)

    return sim.solution
end
