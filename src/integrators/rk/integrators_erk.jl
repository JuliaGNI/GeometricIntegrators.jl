
"Holds the tableau of an explicit Runge-Kutta method."
struct TableauERK{T} <: AbstractTableauERK{T}
    @HeaderTableau

    q::CoefficientsRK{T}

    function TableauERK{T}(q) where {T}
        @assert q.c[1] == 0
        @assert istrilstrict(q.a)
        @assert !(q.s==1 && q.a[1,1] ≠ 0)
        new(q.name, q.o, q.s, q)
    end
end

function TableauERK(q::CoefficientsRK{T}) where {T}
    TableauERK{T}(q)
end

function TableauERK(name::Symbol, order::Int, a::Matrix{T}, b::Vector{T}, c::Vector{T}) where {T}
    TableauERK{T}(CoefficientsRK(name, order, a, b, c))
end

"Read explicit Runge-Kutta tableau from file."
function readTableauERKFromFile(dir::AbstractString, name::AbstractString)
    file = string(dir, "/", name, ".tsv")

#    run(`cat $file`)

    o, s, T = readTableauRKHeaderFromFile(file)

    # TODO Read data in original format (e.g., Rational).
    #      For this we need to save tableaus as jld or hdf5.
#    tab_array = readdlm(file, T)
    tab_array = readdlm(file, comments=true)

    if s == 0
        s = size(tab_array, 1)
    end

    @assert s == size(tab_array, 1) == size(tab_array, 2)-2

    c = tab_array[1:s, 1]
    b = tab_array[1:s, 2]
    a = tab_array[1:s, 3:s+2]

    @info "Reading explicit Runge-Kutta tableau $(name) with $(s) stages and order $(o) from file\n$(file)"

    TableauERK(Symbol(name), o, a, b, c)
end



"Explicit Runge-Kutta integrator."
struct IntegratorERK{DT,TT,FT} <: Integrator{DT,TT}
    equation::ODE{DT,TT,FT}
    tableau::TableauERK{TT}
    Δt::TT

    q::Vector{DT}
    Q::Vector{Vector{DT}}
    V::Vector{Vector{DT}}

    function IntegratorERK{DT,TT,FT}(equation, tableau, Δt) where {DT,TT,FT}
        D = equation.d
        S = tableau.q.s
        Q = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        new(equation, tableau, Δt, zeros(DT,D), Q, V)
    end
end

function IntegratorERK(equation::ODE{DT,TT,FT}, tableau::TableauERK{TT}, Δt::TT) where {DT,TT,FT}
    IntegratorERK{DT,TT,FT}(equation, tableau, Δt)
end

function initialize!(int::IntegratorERK, sol::SolutionODE, m::Int)
    @assert m ≥ 1
    @assert m ≤ sol.ni

    # copy initial conditions from solution
    get_initial_conditions!(sol, int.q[m], m)
end

"Integrate ODE with explicit Runge-Kutta integrator."
function integrate_step!(int::IntegratorERK{DT,TT}, sol::SolutionODE{DT,TT}, m::Int, n::Int) where {DT,TT}
    @assert m ≥ 1
    @assert m ≤ sol.ni

    @assert n ≥ 1
    @assert n ≤ sol.ntime

    local tᵢ::TT
    local yᵢ::DT

    # compute internal stages
    for i in 1:int.tableau.q.s
        @inbounds for k in eachindex(int.Q[i], int.V[i])
            yᵢ = 0
            for j in 1:i-1
                yᵢ += int.tableau.q.a[i,j] * int.V[j][k]
            end
            int.Q[i][k] = int.q[k] + int.Δt * yᵢ
        end
        tᵢ = sol.t[0] + (n-1)*int.Δt + int.Δt * int.tableau.q.c[i]
        int.equation.v(tᵢ, int.Q[i], int.V[i])
    end

    # compute final update
    update_solution!(int.q[m], int.V, int.tableau.q.b, int.Δt)

    # take care of periodic solutions
    cut_periodic_solution!(int.q[m], int.equation.periodicity)

    # copy to solution
    copy_solution!(sol, int.q[m], n, m)
end
