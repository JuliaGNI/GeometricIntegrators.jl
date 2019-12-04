
# general helper functions for integrators.

"""
Create nonlinear solver object for a system of `N` equations with data type `DT`.
The function ``f(x)=0`` to be solved for is determined by a julia function
`function_stages!(x, b, params)`, where `x` is the current solution and `b` is
the output vector, s.th. ``b = f(x)``. `params` are a set of parameters depending
on the equation and integrator that is used.
The solver type is obtained from the config dictionary (`:nls_solver`).
"""
function create_nonlinear_solver(DT, N, params)
    # create solution vector for nonlinear solver
    x = zeros(DT, N)

    # create wrapper function f(x,b) that calls `function_stages!(x, b, params)`
    # with the appropriate `params`
    f = (x,b) -> function_stages!(x, b, params)

    # create nonlinear solver with solver type obtained from config dictionary
    s = get_config(:nls_solver)(x, f)
end


"""
Create a solution vector of type `TwicePrecision{DT}` for a problem with `D` dimensions
and `M` independent initial conditions.
"""
function create_solution_vector(DT, D, M)
    [zeros(TwicePrecision{DT}, D) for i in 1:M]
end


"""
Create a solution vector of type `TwicePrecision{DT}` for a problem with `D` dimensions,
`NS' sample paths, and `NI` independent initial conditions.
"""
function create_solution_vector(DT, D, NS, NI)
    [zeros(TwicePrecision{DT}, D) for i in 1:NS, j in 1:NI]
end


"""
Create a vector of S solution vectors of type DT to store the solution of S
internal stages for a problem with `D` dimensions.
"""
function create_internal_stage_vector(DT, D, S)
    [zeros(DT,D) for i in 1:S]
end


"""
Create a vector of S solution matrices of type DT to store the solution of S
internal stages for a problem with `DxM` dimensions.
"""
function create_internal_stage_vector(DT, D, M, S)
    [zeros(DT,D,M) for i in 1:S]
end


"""
Create a vector of S+1 solution vectors of type DT to store the solution of S
internal stages and the solution of the previous timestep for a problem with `D`
    dimensions.
"""
function create_internal_stage_vector_with_zero(DT, D, S)
    a = OffsetArray{Vector{DT}}(undef, 0:S)

    for i in 0:S
        a[i] = zeros(DT,D)
    end

    return a
end


function check_solution_dimension_asserts(sol::Solution, m::Int, n::Int)
    @assert m ≥ 1
    @assert m ≤ sol.ni

    @assert n ≥ 1
    @assert n ≤ sol.ntime
end


function update_solution!(x::Vector{T}, xₑᵣᵣ::Vector{T}, ẋ::Matrix{T}, b::Vector{T}, Δt::T) where {T}
    @assert length(x) == length(xₑᵣᵣ)
    @assert length(x) == size(ẋ, 1)
    @assert length(b) == size(ẋ, 2)

    for k in axes(ẋ, 1)
        for i in axes(ẋ, 2)
            x[k], xₑᵣᵣ[k] = compensated_summation(Δt * b[i] * ẋ[k,i], x[k], xₑᵣᵣ[k])
        end
    end
end

function update_solution!(x::Vector{T}, xₑᵣᵣ::Vector{T}, ẋ::Vector{Vector{T}}, b::Vector{T}, Δt::T) where {T}
    @assert length(b) == length(ẋ)
    @assert length(x) == length(ẋ[1])
    @assert length(x) == length(xₑᵣᵣ)

    for i in eachindex(ẋ)
        for k in eachindex(ẋ[i])
            x[k], xₑᵣᵣ[k] = compensated_summation(Δt * b[i] * ẋ[i][k], x[k], xₑᵣᵣ[k])
        end
    end
end

function update_solution!(x::Union{Vector{T}, Vector{TwicePrecision{T}}}, ẋ::Matrix{T}, b::Vector{T}, Δt::T) where {T}
    @assert length(x) == size(ẋ, 1)
    @assert length(b) == size(ẋ, 2)

    local Δx::eltype(x)

    for k in axes(ẋ, 1)
        Δx = 0
        for i in axes(ẋ, 2)
            Δx += b[i] * ẋ[k,i]
        end
        x[k] += Δt * Δx
    end
end


# For stochastic Runge-Kutta methods (SIRK and WIRK)
# x - the solution vector to be updated
# Vx - the matrix containing the drift vector evaluated at the internal stages v(Q_i) (SIRK) or v(Q0_i) (WIRK)
# Bx - the array containing the diffusion matrix evaluated at the internal stages B(Q_i) (SIRK) or B(Q1^(l)_i) (WIRK)
# bdrift - the Runge-Kutta coefficients for the drift part
# bdiff - the Runge-Kutta coefficients for the diffusion part
# Δt - the time step
# ΔW - the increments of the Brownian motion (SFIRK) or the increments represented by the random variables \hat I^(k) (WFIRK)
function update_solution!(x::Union{Vector{T}, Vector{TwicePrecision{T}}}, Vx::Vector{Vector{T}}, Bx::Vector{Matrix{T}}, bdrift::Vector{T}, bdiff::Vector{T}, Δt::T, ΔW::Vector{T}) where {T}

    @assert length(bdrift) == length(bdiff) == length(Vx) == length(Bx)

    for i in eachindex(Vx, Bx)
        @assert length(x) == length(Vx[i]) == size(Bx[i], 1)
        @assert length(ΔW)== size(Bx[i], 2)
    end

    local Δx::eltype(x)

    # Contribution from the drift part
    for k in eachindex(x)
        Δx = 0.
        for i in eachindex(Vx)
            Δx += bdrift[i] * Vx[i][k]
        end
        x[k] += Δt * Δx
    end

    local Δy = zero(ΔW)

    # Contribution from the diffusion part
    for k in eachindex(x)
        Δy .= 0
        for i in eachindex(Bx)
            Δy .+= bdiff[i] .* Bx[i][k,:]
        end
        x[k] += dot(Δy,ΔW)
    end
end


# For stochastic Runge-Kutta methods
# x - the solution vector to be updated
# Vx - the matrix containing the drift vector evaluated at the internal stages v(Q_i)
# Bx - the array containing the diffusion matrix evaluated at the internal stages B(Q_i)
# bdrift - the Runge-Kutta coefficients for the drift part
# bdiff - the Runge-Kutta coefficients for the ΔW terms of the diffusion part
# bdiff2- the Runge-Kutta coefficients for the ΔZ terms of the diffusion part
# Δt - the time step
# ΔW - the increments of the Brownian motion
# ΔZ - the integrals of the increments of the Brownian motion
function update_solution!(x::Union{Vector{T}, Vector{TwicePrecision{T}}}, Vx::Vector{Vector{T}}, Bx::Vector{Matrix{T}}, bdrift::Vector{T}, bdiff::Vector{T}, bdiff2::Vector{T}, Δt::T, ΔW::Vector{T}, ΔZ::Vector{T}) where {T}

    @assert length(bdrift) == length(bdiff) == length(bdiff2) == length(Vx) == length(Bx)

    for i in eachindex(Vx, Bx)
        @assert length(x) == length(Vx[i]) == size(Bx[i], 1)
        @assert length(ΔW) == length(ΔZ) == size(Bx[i], 2)
    end

    local Δx::eltype(x)

    # Contribution from the drift part
    for k in eachindex(x)
        Δx = 0.
        for i in eachindex(Vx)
            Δx += bdrift[i] * Vx[i][k]
        end
        x[k] += Δt * Δx
    end

    local Δy = zero(ΔW)

    # Contribution from the diffusion part (ΔW terms)
    for k in eachindex(x)
        Δy .= 0
        for i in eachindex(Bx)
            Δy .+= bdiff[i] .* Bx[i][k,:]
        end
        x[k] += dot(Δy,ΔW)
    end

    # Contribution from the diffusion part (ΔZ terms)
    for k in eachindex(x)
        Δy .= 0
        for i in eachindex(Bx)
            Δy .+= bdiff2[i] .* Bx[i][k,:]
        end
        x[k] += dot(Δy,ΔZ)/Δt
    end
end


# For stochastic partitioned Runge-Kutta methods
# q, p - the solution vector to be updated
# Vqp, Fqp - the matrix containing the drift vectors evaluated at the internal stages v(Q_i), f(Q_i)
# Bqp, Gqp - the array containing the diffusion matrices evaluated at the internal stages B(Q_i), G(Q_i)
# bqdrift, bpdrift - the Runge-Kutta coefficients for the drift parts of the q and p equations
# bqdiff, bpdiff - the Runge-Kutta coefficients for the diffusion parts of the q and p equations
# Δt - the time step
# ΔW - the increments of the Brownian motion
function update_solution!(q::Vector{VT}, p::Vector{VT},
                          Vqp::Vector{Vector{T}}, Fqp::Vector{Vector{T}},
                          Bqp::Vector{Matrix{T}}, Gqp::Vector{Matrix{T}},
                          bqdrift::Vector{T}, bqdiff::Vector{T},
                          bpdrift::Vector{T}, bpdiff::Vector{T},
                          Δt::T, ΔW::Vector{T}) where {T, VT <: Union{T, TwicePrecision{T}}}

    @assert length(bqdrift) == length(bqdiff) == length(bpdrift) == length(bpdiff) == length(Vqp) == length(Fqp) == length(Bqp) == length(Gqp)

    for i in eachindex(Vqp, Fqp, Bqp, Gqp)
        @assert length(q) == length(p) == length(Vqp[i]) == length(Fqp[i]) == size(Bqp[i], 1) == size(Gqp[i], 1)
        @assert length(ΔW) == size(Bqp[i], 2) == size(Gqp[i], 2)
    end

    local Δq::eltype(q)
    local Δp::eltype(p)

    # Contribution from the drift part
    for k in eachindex(q, p)
        Δq = 0.
        Δp = 0.
        for i in eachindex(Vqp, Fqp)
            Δq += bqdrift[i] * Vqp[i][k]
            Δp += bpdrift[i] * Fqp[i][k]
        end
        q[k] += Δt * Δq
        p[k] += Δt * Δp
    end

    local Δy = zero(ΔW)
    local Δz = zero(ΔW)

    # Contribution from the diffusion part
    for k in eachindex(q, p)
        Δy .= 0
        Δz .= 0

        for i in eachindex(Bqp, Gqp)
            Δy .+= bqdiff[i] .* Bqp[i][k,:]
            Δz .+= bpdiff[i] .* Gqp[i][k,:]
        end

        q[k] += dot(Δy,ΔW)
        p[k] += dot(Δz,ΔW)
    end
end


# For stochastic split partitioned Runge-Kutta methods
# q, p - the solution vector to be updated
# Vqp, Fqp1, Fqp2 - the matrix containing the drift vectors evaluated at the internal stages v(Q_i), fi(Q_i)
# Bqp, Gqp1, Gqp2 - the array containing the diffusion matrices evaluated at the internal stages B(Q_i), Gi(Q_i)
# bqdrift, bpdrift1, bpdrift2 - the Runge-Kutta coefficients for the drift parts of the q and p equations
# bqdiff, bpdiff1, bpdiff2 - the Runge-Kutta coefficients for the diffusion parts of the q and p equations
# Δt - the time step
# ΔW - the increments of the Brownian motion
function update_solution!(q::Union{Vector{T}, Vector{TwicePrecision{T}}}, p::Union{Vector{T}, Vector{TwicePrecision{T}}},
                            Vqp::Vector{Vector{T}}, Fqp1::Vector{Vector{T}}, Fqp2::Vector{Vector{T}},
                            Bqp::Vector{Matrix{T}}, Gqp1::Vector{Matrix{T}}, Gqp2::Vector{Matrix{T}},
                            bqdrift::Vector{T}, bqdiff::Vector{T},
                            bpdrift1::Vector{T}, bpdrift2::Vector{T},
                            bpdiff1::Vector{T}, bpdiff2::Vector{T},
                            Δt::T, ΔW::Vector{T}) where {T}

   @assert length(bqdrift) == length(bqdiff) == length(bpdrift1) == length(bpdrift2) == length(bpdiff1) == length(bpdiff2) == length(Vqp) == length(Fqp1) == length(Fqp2) == length(Bqp) == length(Gqp1) == length(Gqp2)

   for i in eachindex(Vqp, Fqp1, Fqp2, Bqp, Gqp1, Gqp2)
       @assert length(q) == length(p) == length(Vqp[i]) == length(Fqp1[i]) == length(Fqp2[i]) == size(Bqp[i], 1) == size(Gqp1[i], 1) == size(Gqp2[i], 1)
       @assert length(ΔW) == size(Bqp[i], 2) == size(Gqp1[i], 2) == size(Gqp2[i], 2)
   end

   local Δq::eltype(q)
   local Δp::eltype(p)

   # Contribution from the drift part
   for k in eachindex(q, p)
       Δq = 0.
       Δp = 0.
       for i in eachindex(Vqp, Fqp1, Fqp2)
           Δq += bqdrift[i] * Vqp[i][k]
           Δp += bpdrift1[i] * Fqp1[i][k] + bpdrift2[i] * Fqp2[i][k]
       end
       q[k] += Δt * Δq
       p[k] += Δt * Δp
   end

   local Δy = zero(ΔW)
   local Δz = zero(ΔW)

   # Contribution from the diffusion part
   for k in eachindex(q, p)
       Δy .= 0
       Δz .= 0

       for i in eachindex(Bqp, Gqp1, Gqp2)
           Δy .+= bqdiff[i]  .* Bqp[i][k,:]
           Δz .+= bpdiff1[i] .* Gqp1[i][k,:] .+ bpdiff2[i] .* Gqp2[i][k,:]
       end

       q[k] += dot(Δy,ΔW)
       p[k] += dot(Δz,ΔW)
   end
end


# For weak Runge-Kutta methods WERK
# x - the solution vector to be updated
# Vx - the matrix containing the drift vector evaluated at the internal stages v(Q_i)
# Bx1 - the array containing the diffusion matrix evaluated at the internal stages H^(l)_i, such that Bx1[:,l,i] is evaluated at H^(l)_i
# Bx2 - the array containing the diffusion matrix evaluated at the internal stages \hat H^(l)_i, such that Bx2[:,l,i] is evaluated at \hat H^(l)_i
# alpha - the Runge-Kutta coefficients for the drift part
# beta1 - the Runge-Kutta coefficients for the diffusion term with the random increments
# beta2 - the Runge-Kutta coefficients for the second diffusion term
# Δt - the time step
# ΔW - the increments of the Brownian motion represented by the random variables \hat I^(k)
function update_solution!(x::Union{Vector{T}, Vector{TwicePrecision{T}}}, Vx::Vector{Vector{T}}, Bx1::Vector{Matrix{T}}, Bx2::Vector{Matrix{T}}, alpha::Vector{T}, beta1::Vector{T}, beta2::Vector{T}, Δt::T, ΔW::Vector{T}) where {T}
    @assert length(alpha) == length(beta1) == length(beta2) == length(Vx) == length(Bx1) == length(Bx2)

    for i in eachindex(Vx, Bx1, Bx2)
        @assert length(x) == length(Vx[i]) == size(Bx1[i], 1) == size(Bx2[i], 1)
        @assert length(ΔW)== size(Bx1[i], 2) == size(Bx2[i], 2)
    end

    local Δx::eltype(x)

    # Contribution from the drift part
    for k in eachindex(x)
        Δx = 0.
        for i in eachindex(Vx)
            Δx += alpha[i] * Vx[i][k]
        end
        x[k] += Δt * Δx
    end

    local Δy = zero(ΔW)

    # Contribution from the diffusion term with the random variables I^(k)_i
    for k in eachindex(x)
        Δy .= 0
        for i in eachindex(Bx1)
            Δy .+= beta1[i] .* Bx1[i][k,:]
        end
        x[k] += dot(Δy,ΔW)
    end

    # Contribution from the second diffusion term
    for k in eachindex(x)
        Δx= 0.
        for i in eachindex(Bx2)
            for l in axes(Bx2[i], 2)
                Δx += beta2[i] * Bx2[i][k,l]
            end
        end
        x[k] += sqrt(Δt)*Δx
    end
end


function update_solution!(x::Union{Vector{T}, Vector{TwicePrecision{T}}}, ẋ::Vector{Vector{T}}, b::Vector{T}, Δt::T) where {T}
    @assert length(b) == length(ẋ)
    @assert length(x) == length(ẋ[1])

    for i in eachindex(ẋ)
        for k in eachindex(ẋ[i])
            x[k] += Δt * b[i] * ẋ[i][k]
        end
    end
end

function update_solution!(x::Vector{T}, xₑᵣᵣ::Vector{T}, ẋ::Union{Matrix{T},Vector{Vector{T}}}, b::Vector{T}, b̂::Vector, Δt::T) where {T}
    update_solution!(x, xₑᵣᵣ, ẋ, b, Δt)
    update_solution!(x, xₑᵣᵣ, ẋ, b̂, Δt)
end

function update_solution!(x::Union{Vector{T}, Vector{TwicePrecision{T}}}, ẋ::Union{Matrix{T},Vector{Vector{T}}}, b::Vector{T}, b̂::Vector, Δt::T) where {T}
    update_solution!(x, ẋ, b, Δt)
    update_solution!(x, ẋ, b̂, Δt)
end

# For stochastic Runge-Kutta methods
function update_solution!(x::Union{Vector{T}, Vector{TwicePrecision{T}}}, Vx::Vector{Vector{T}}, Bx::Vector{Matrix{T}},
                            bdrift::Vector{T}, b̂drift::Vector, bdiff::Vector{T}, b̂diff::Vector, Δt::T, ΔW::Vector{T}) where {T}
    update_solution!(x, Vx, Bx, bdrift, bdiff, Δt, ΔW)
    update_solution!(x, Vx, Bx, b̂drift, b̂diff, Δt, ΔW)
end

# For stochastic partitioned Runge-Kutta methods
function update_solution!(q::Union{Vector{T}, Vector{TwicePrecision{T}}}, p::Union{Vector{T}, Vector{TwicePrecision{T}}},
                            Vqp::Vector{Vector{T}}, Fqp::Vector{Vector{T}}, Bqp::Vector{Matrix{T}}, Gqp::Vector{Matrix{T}},
                            bqdrift::Vector{T}, b̂qdrift::Vector, bqdiff::Vector{T}, b̂qdiff::Vector,
                            bpdrift::Vector{T}, b̂pdrift::Vector, bpdiff::Vector{T}, b̂pdiff::Vector, Δt::T, ΔW::Vector{T}) where {T}
    update_solution!(q, p, Vqp, Fqp, Bqp, Gqp, bqdrift, bqdiff, bpdrift, bpdiff, Δt, ΔW)
    update_solution!(q, p, Vqp, Fqp, Bqp, Gqp, b̂qdrift, b̂qdiff, b̂pdrift, b̂pdiff, Δt, ΔW)
end

# For stochastic split partitioned Runge-Kutta methods
function update_solution!(q::Union{Vector{T}, Vector{TwicePrecision{T}}}, p::Union{Vector{T}, Vector{TwicePrecision{T}}},
                            Vqp::Vector{Vector{T}}, Fqp1::Vector{Vector{T}}, Fqp2::Vector{Vector{T}},
                            Bqp::Vector{Matrix{T}}, Gqp1::Vector{Matrix{T}}, Gqp2::Vector{Matrix{T}},
                            bqdrift::Vector{T}, b̂qdrift::Vector, bqdiff::Vector{T}, b̂qdiff::Vector,
                            bpdrift1::Vector{T}, b̂pdrift1::Vector, bpdrift2::Vector{T}, b̂pdrift2::Vector,
                            bpdiff1::Vector{T}, b̂pdiff1::Vector, bpdiff2::Vector{T}, b̂pdiff2::Vector,
                            Δt::T, ΔW::Vector{T}) where {T}
    update_solution!(q, p, Vqp, Fqp1, Fqp2, Bqp, Gqp1, Gqp2, bqdrift, bqdiff, bpdrift1, bpdrift2, bpdiff1, bpdiff2, Δt, ΔW)
    update_solution!(q, p, Vqp, Fqp1, Fqp2, Bqp, Gqp1, Gqp2, b̂qdrift, b̂qdiff, b̂pdrift1, b̂pdrift2, b̂pdiff1, b̂pdiff2, Δt, ΔW)
end

 # For weak Runge-Kutta methods WERK
 function update_solution!(x::Union{Vector{T}, Vector{TwicePrecision{T}}}, Vx::Vector{Vector{T}}, Bx1::Vector{Matrix{T}}, Bx2::Vector{Matrix{T}},
                             alpha::Vector{T}, âlpha::Vector, beta1::Vector{T}, b̂eta1::Vector, beta2::Vector{T}, b̂eta2::Vector, Δt::T, ΔW::Vector{T}) where {T}
     update_solution!(x, Vx, Bx1, Bx2, alpha, beta1, beta2, Δt, ΔW)
     update_solution!(x, Vx, Bx1, Bx2, âlpha, b̂eta1, b̂eta2, Δt, ΔW)
 end


function update_multiplier!(λ::SolutionVector{T}, Λ::Vector{Vector{T}}, b::Vector{T}) where {T}
    @assert length(λ) == size(Λ, 1)
    local t::T
    @inbounds for i in eachindex(λ)
        t = zero(T)
        for j=eachindex(b)
            t += b[j] * Λ[j][i]
        end
        λ[i] = t
    end
    nothing
end


function CommonFunctions.cut_periodic_solution!(x::Vector{T}, periodicity::Vector{T}) where {T}
    @assert length(x) == length(periodicity)

    for k in eachindex(x, periodicity)
        if periodicity[k] ≠ 0
            while x[k] < 0
                x[k] += periodicity[k]
            end
            while x[k] ≥ periodicity[k]
                x[k] -= periodicity[k]
            end
        end
    end
end

function CommonFunctions.cut_periodic_solution!(x::Vector{TwicePrecision{T}}, periodicity::Vector{T}) where {T}
    @assert length(x) == length(periodicity)

    for k in eachindex(x, periodicity)
        if periodicity[k] ≠ 0
            while x[k].hi < 0
                x[k] += periodicity[k]
            end
            while x[k].hi ≥ periodicity[k]
                x[k] -= periodicity[k]
            end
        end
    end
end

function CommonFunctions.cut_periodic_solution!(x::Vector{T}, periodicity::Vector{T}, shift::Vector{T}) where {T}
    @assert length(x) == length(periodicity)
    shift .= 0
    for k in eachindex(x, periodicity, shift)
        if periodicity[k] ≠ 0
            while x[k] + shift[k] < 0
                shift[k] += periodicity[k]
            end
            while x[k] + shift[k] ≥ periodicity[k]
                shift[k] -= periodicity[k]
            end
        end
    end
end

function CommonFunctions.cut_periodic_solution!(x::Vector{TwicePrecision{T}}, periodicity::Vector{T}, shift::Vector{T}) where {T}
    @assert length(x) == length(periodicity)
    shift .= 0
    for k in eachindex(x, periodicity, shift)
        if periodicity[k] ≠ 0
            while x[k].hi + shift[k] < 0
                shift[k] += periodicity[k]
            end
            while x[k].hi + shift[k] ≥ periodicity[k]
                shift[k] -= periodicity[k]
            end
        end
    end
end
