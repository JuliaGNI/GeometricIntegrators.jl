
abstract type SolutionPDAE{dType, tType, N} <: DeterministicSolution{dType, tType, N} end


"Serial Solution of a partitioned differential algebraic equation."
mutable struct SSolutionPDAE{dType, tType, N} <: SolutionPDAE{dType, tType, N}
    nd::Int
    nm::Int
    nt::Int
    ni::Int
    t::TimeSeries{tType}
    q::SDataSeries{dType,N}
    p::SDataSeries{dType,N}
    λ::SDataSeries{dType,N}
    ntime::Int
    nsave::Int
    counter::Int
end

function SSolutionPDAE(equation::Union{IODE{DT,TT},VODE{DT,TT},PDAE{DT,TT},IDAE{DT,TT}}, Δt::TT, ntime::Int, nsave::Int=1) where {DT,TT}
    N  = equation.n > 1 ? 3 : 2
    nd = equation.d
    nm = equation.m
    ni = equation.n
    nt = div(ntime, nsave)

    @assert DT <: Number
    @assert TT <: Real
    @assert nd > 0
    @assert nm > 0
    @assert ni > 0
    @assert nsave > 0
    @assert ntime == 0 || ntime ≥ nsave
    @assert mod(ntime, nsave) == 0

    t = TimeSeries{TT}(nt, Δt, nsave)
    q = SDataSeries(DT, nd, nt, ni)
    p = SDataSeries(DT, nd, nt, ni)
    λ = SDataSeries(DT, nm, nt, ni)
    s = SSolutionPDAE{DT,TT,N}(nd, nm, nt, ni, t, q, p, λ, ntime, nsave, 0)
    set_initial_conditions!(s, equation)
    return s
end

function SSolutionPDAE(t::TimeSeries{TT}, q::SDataSeries{DT,N}, p::SDataSeries{DT,N}, λ::SDataSeries{DT,N}, ntime::Int, nsave::Int) where {DT,TT,N}
    @assert q.nd == p.nd
    @assert q.nt == p.nt == λ.nt
    @assert q.ni == p.ni == λ.ni

    # extract parameters
    nd = q.nd
    nm = λ.nd
    ni = q.ni
    nt = t.n

    # create solution
    SSolutionPDAE{DT,TT,N}(nd, nm, nt, ni, t, q, p, λ, ntime, nsave, 0)
end

function SSolutionPDAE(file::String)
    # open HDF5 file
    info("Reading HDF5 file ", file)
    h5 = h5open(file, "r")

    # read attributes
    ntime = read(attrs(h5)["ntime"])
    nsave = read(attrs(h5)["nsave"])

    # reading data arrays
    t = TimeSeries(read(h5["t"]), nsave)
    q = SDataSeries(read(h5["q"]))
    p = SDataSeries(read(h5["p"]))
    λ = SDataSeries(read(h5["λ"]))

    # create solution
    SSolutionPDAE(t, q, p, λ, ntime, nsave)
end


"Parallel Solution of a partitioned differential algebraic equation."
mutable struct PSolutionPDAE{dType, tType, N} <: SolutionPDAE{dType, tType, N}
    nd::Int
    nm::Int
    nt::Int
    ni::Int
    t::TimeSeries{tType}
    q::PDataSeries{dType,N}
    p::PDataSeries{dType,N}
    λ::PDataSeries{dType,N}
    ntime::Int
    nsave::Int
    counter::Int
end

function PSolutionPDAE(equation::Union{IODE{DT,TT},VODE{DT,TT},PDAE{DT,TT},IDAE{DT,TT}}, Δt::TT, ntime::Int, nsave::Int=1) where {DT,TT}
    N  = equation.n > 1 ? 3 : 2
    nd = equation.d
    nm = equation.m
    ni = equation.n
    nt = div(ntime, nsave)

    @assert DT <: Number
    @assert TT <: Real
    @assert nd > 0
    @assert nm > 0
    @assert ni > 0
    @assert nsave > 0
    @assert ntime == 0 || ntime ≥ nsave
    @assert mod(ntime, nsave) == 0

    t = TimeSeries{TT}(nt, Δt, nsave)
    q = PDataSeries(DT, nd, nt, ni)
    p = PDataSeries(DT, nd, nt, ni)
    λ = PDataSeries(DT, nm, nt, ni)
    s = PSolutionPDAE{DT,TT,N}(nd, nm, nt, ni, t, q, p, λ, ntime, nsave, 0)
    set_initial_conditions!(s, equation)
    return s
end

function PSolutionPDAE(t::TimeSeries{TT}, q::PDataSeries{DT,N}, p::PDataSeries{DT,N}, λ::PDataSeries{DT,N}, ntime::Int, nsave::Int) where {DT,TT,N}
    @assert q.nd == p.nd
    @assert q.nt == p.nt == λ.nt
    @assert q.ni == p.ni == λ.ni

    # extract parameters
    nd = q.nd
    nm = λ.nd
    ni = q.ni
    nt = t.n

    # create solution
    PSolutionPDAE{DT,TT,N}(nd, nm, nt, ni, t, q, p, λ, ntime, nsave, 0)
end

function PSolutionPDAE(file::String)
    # open HDF5 file
    info("Reading HDF5 file ", file)
    h5 = h5open(file, "r")

    # read attributes
    ntime = read(attrs(h5)["ntime"])
    nsave = read(attrs(h5)["nsave"])

    # reading data arrays
    t = TimeSeries(read(h5["t"]), nsave)
    q = PDataSeries(read(h5["q"]))
    p = PDataSeries(read(h5["p"]))
    λ = PDataSeries(read(h5["λ"]))

    # create solution
    PSolutionPDAE(t, q, p, λ, ntime, nsave)
end


time(sol::SolutionPDAE)  = sol.t.t
ntime(sol::SolutionPDAE) = sol.ntime
nsave(sol::SolutionPDAE) = sol.nsave


function set_initial_conditions!(sol::SolutionPDAE{DT,TT}, equ::Union{IODE{DT,TT},VODE{DT,TT},PDAE{DT,TT},IDAE{DT,TT}}) where {DT,TT}
    set_initial_conditions!(sol, equ.t₀, equ.q₀, equ.p₀, equ.λ₀)
end

function set_initial_conditions!(sol::SolutionPDAE{DT,TT}, t₀::TT, q₀::Union{Array{DT}, Array{TwicePrecision{DT}}}, p₀::Union{Array{DT}, Array{TwicePrecision{DT}}}, λ₀::Union{Array{DT}, Array{TwicePrecision{DT}}}) where {DT,TT}
    set_data!(sol.q, q₀, 0)
    set_data!(sol.p, p₀, 0)
    set_data!(sol.λ, λ₀, 0)
    compute_timeseries!(sol.t, t₀)
end

function get_initial_conditions!(sol::SolutionPDAE{DT,TT}, q::Union{Vector{DT}, Vector{TwicePrecision{DT}}}, p::Union{Vector{DT}, Vector{TwicePrecision{DT}}}, λ::Union{Vector{DT}, Vector{TwicePrecision{DT}}}, k) where {DT,TT}
    get_data!(sol.q, q, 0, k)
    get_data!(sol.p, p, 0, k)
    get_data!(sol.λ, λ, 0, k)
end

function get_initial_conditions!(sol::SolutionPDAE{DT,TT}, q::Union{Vector{DT}, Vector{TwicePrecision{DT}}}, p::Union{Vector{DT}, Vector{TwicePrecision{DT}}}, k) where {DT,TT}
    get_data!(sol.q, q, 0, k)
    get_data!(sol.p, p, 0, k)
end

function copy_solution!(sol::SolutionPDAE{DT,TT}, q::Union{Vector{DT}, Vector{TwicePrecision{DT}}}, p::Union{Vector{DT}, Vector{TwicePrecision{DT}}}, λ::Union{Vector{DT}, Vector{TwicePrecision{DT}}}, n, k) where {DT,TT}
    if mod(n, sol.nsave) == 0
        j = div(n, sol.nsave)
        set_data!(sol.q, q, j, k)
        set_data!(sol.p, p, j, k)
        set_data!(sol.λ, λ, j, k)
        sol.counter += 1
    end
end

function copy_solution!(sol::SolutionPDAE{DT,TT}, q::Union{Vector{DT}, Vector{TwicePrecision{DT}}}, p::Union{Vector{DT}, Vector{TwicePrecision{DT}}}, n, k) where {DT,TT}
    if mod(n, sol.nsave) == 0
        j = div(n, sol.nsave)
        set_data!(sol.q, q, j, k)
        set_data!(sol.p, p, j, k)
        sol.counter += 1
    end
end

function reset!(sol::SolutionPDAE)
    reset!(sol.q)
    reset!(sol.p)
    reset!(sol.λ)
    compute_timeseries!(sol.t, sol.t[end])
    sol.counter = 0
end


"Creates HDF5 file and initialises datasets for PDAE solution object with single initial condition."
function create_hdf5(solution::SolutionPDAE{DT,TT,2}, file::AbstractString) where {DT,TT}
    # create HDF5 file and save attributes and common parameters
    h5 = createHDF5(solution, file)

    # create datasets
    q = d_create(h5, "q", datatype(Float64), dataspace(solution.nd, solution.nt+1), "chunk", (solution.nd,1))
    p = d_create(h5, "p", datatype(Float64), dataspace(solution.nd, solution.nt+1), "chunk", (solution.nd,1))
    λ = d_create(h5, "λ", datatype(Float64), dataspace(solution.nm, solution.nt+1), "chunk", (solution.nm,1))

    # copy initial conditions
    q[1:solution.nd, 1] = Array{Float64}(solution.q.d[1:solution.nd, 1])
    p[1:solution.nd, 1] = Array{Float64}(solution.p.d[1:solution.nd, 1])
    λ[1:solution.nm, 1] = Array{Float64}(solution.λ.d[1:solution.nm, 1])

    return h5
end

"Creates HDF5 file and initialises datasets for PDAE solution object with multiple initial conditions."
function create_hdf5(solution::SolutionPDAE{DT,TT,3}, file::AbstractString) where {DT,TT}
    # create HDF5 file and save attributes and common parameters
    h5 = createHDF5(solution, file)

    # create datasets
    q = d_create(h5, "q", datatype(Float64), dataspace(solution.nd, solution.nt+1, solution.ni), "chunk", (solution.nd,1,1))
    p = d_create(h5, "p", datatype(Float64), dataspace(solution.nd, solution.nt+1, solution.ni), "chunk", (solution.nd,1,1))
    λ = d_create(h5, "λ", datatype(Float64), dataspace(solution.nm, solution.nt+1, solution.ni), "chunk", (solution.nm,1,1))

    # copy initial conditions
    q[1:solution.nd, 1, 1:solution.ni] = Array{Float64}(solution.q.d[1:solution.nd, 1, 1:solution.ni])
    p[1:solution.nd, 1, 1:solution.ni] = Array{Float64}(solution.p.d[1:solution.nd, 1, 1:solution.ni])
    λ[1:solution.nm, 1, 1:solution.ni] = Array{Float64}(solution.λ.d[1:solution.nm, 1, 1:solution.ni])

    return h5
end

"Append solution to HDF5 file."
function CommonFunctions.write_to_hdf5(solution::SolutionPDAE{DT,TT,2}, h5::HDF5.HDF5File, offset=0) where {DT,TT}
    # set convenience variables and compute ranges
    d  = solution.nd
    m  = solution.nm
    n  = solution.nt
    j1 = offset+2
    j2 = offset+1+n

    # copy data from solution to HDF5 dataset
    h5["q"][1:d, j1:j2] = Array{Float64}(solution.q.d[1:d, 2:n+1])
    h5["p"][1:d, j1:j2] = Array{Float64}(solution.p.d[1:d, 2:n+1])
    h5["λ"][1:m, j1:j2] = Array{Float64}(solution.λ.d[1:m, 2:n+1])

    return nothing
end

"Append solution to HDF5 file."
function CommonFunctions.write_to_hdf5(solution::SolutionPDAE{DT,TT,3}, h5::HDF5.HDF5File, offset=0) where {DT,TT}
    # set convenience variables and compute ranges
    d  = solution.nd
    m  = solution.nm
    n  = solution.nt
    i  = solution.ni
    j1 = offset+2
    j2 = offset+1+n

    # copy data from solution to HDF5 dataset
    h5["q"][1:d, j1:j2, 1:i] = Array{Float64}(solution.q.d[1:d, 2:n+1, 1:i])
    h5["p"][1:d, j1:j2, 1:i] = Array{Float64}(solution.p.d[1:d, 2:n+1, 1:i])
    h5["λ"][1:m, j1:j2, 1:i] = Array{Float64}(solution.λ.d[1:m, 2:n+1, 1:i])

    return nothing
end
