@doc raw"""
Splitting integrator for the solution of initial value problems
```math
\dot{q} (t) = v(t, q(t)) , \qquad q(t_{0}) = q_{0} ,
```
whose vector field ``v`` is given as a sum of vector fields
```math
v (t) = v_1 (t) + ... + v_r (t) .
```

`IntegratorSplitting` has two constructors:
```julia
IntegratorSplitting{DT,D}(solutions::Tuple, f::Vector{Int}, c::Vector, Δt)
IntegratorSplitting(equation::SODE, tableau::AbstractTableauSplitting, Δt)
```
In the first constructor, `DT` is the data type of the state vector and `D`
the dimension of the system. In the second constructor, this information
is extracted from the equation. 
The tuple `solutions` contains functions implementing the flow (exact solution)
of the vector fields `v_i`. The vectors `f` and `c` define the actual splitting
method: `f` is a vector of indices of the flows in the split equation to be
solved and `c` is a vector of the same size `f` that contains the coefficients
for each splitting step, i.e., the resulting integrator has the form
```math
\varphi_{\tau} = \phi_{c[s] \tau}^{v_{f[s]}} \circ \dotsc \circ \phi_{c[2] \tau}^{v_{f[2]}} \circ \phi_{c[1] \tau}^{v_{f[1]}} .
```
In the second constructor, these vectors are constructed from the tableau and
the equation.

"""
struct IntegratorSplitting{DT, TT, D, S, QT <: Tuple} <: ODEIntegrator{DT,TT}
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

function IntegratorSplitting(problem::SODEProblem{DT}, tableau::ST) where {DT, TT, ST <: AbstractTableauSplitting{TT}}
    @assert hassolution(problem)
    IntegratorSplitting{DT, ndims(problem)}(solutions(problem).q, get_splitting_coefficients(length(equation(problem).q), tableau)..., timestep(problem))
end


@inline Base.ndims(::IntegratorSplitting{DT,TT,D}) where {DT,TT,D} = D
@inline GeometricBase.timestep(int::IntegratorSplitting) = int.Δt


function integrate_step!(int::IntegratorSplitting{DT,TT}, sol::AtomicSolutionODE{DT,TT}) where {DT,TT}
    local cᵢ::TT
    local tᵢ::TT

    # compute splitting steps
    for i in eachindex(int.f, int.c)
        if int.c[i] ≠ zero(TT)
            cᵢ = timestep(int) * int.c[i]
            tᵢ = sol.t̄ + cᵢ

            # reset atomic solution
            reset!(sol)

            # compute new solution
            int.q[int.f[i]](tᵢ, sol.q̄, sol.q, cᵢ)
        end
    end
end
