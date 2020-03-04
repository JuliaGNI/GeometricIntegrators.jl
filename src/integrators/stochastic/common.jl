# For stochastic Runge-Kutta methods (SIRK and WIRK)
# x - the solution vector to be updated
# V - the matrix containing the drift vector evaluated at the internal stages v(Q_i) (SIRK) or v(Q0_i) (WIRK)
# B - the array containing the diffusion matrix evaluated at the internal stages B(Q_i) (SIRK) or B(Q1^(l)_i) (WIRK)
# bdrift - the Runge-Kutta coefficients for the drift part
# bdiff  - the Runge-Kutta coefficients for the diffusion part
# Δt - the time step
# ΔW - the increments of the Brownian motion (SFIRK) or the increments represented by the random variables \hat I^(k) (WFIRK)
function update_solution!(x::SolutionVector{T}, V::Vector{Vector{T}}, B::Vector{Matrix{T}},
                            bdrift::Vector{T}, bdiff::Vector{T}, Δt::T, ΔW::Vector{T}, Δy::Vector{T}=zero(ΔW)) where {T}

    @assert length(bdrift) == length(bdiff) == length(V) == length(B)

    for i in eachindex(V, B)
        @assert length(x) == length(V[i]) == size(B[i], 1)
        @assert length(ΔW)== size(B[i], 2)
    end

    local Δx::T

    # Contribution from the drift part
    for k in eachindex(x)
        Δx = 0
        for i in eachindex(bdrift, V)
            Δx += bdrift[i] * V[i][k]
        end
        x[k] += Δt * Δx
    end

    # Contribution from the diffusion part
    for k in eachindex(x)
        Δy .= 0
        for i in eachindex(bdiff, B)
            for l in eachindex(Δy)
                Δy[l] += bdiff[i] * B[i][k,l]
            end
        end
        x[k] += dot(Δy,ΔW)
    end
end


# For stochastic Runge-Kutta methods
# x - the solution vector to be updated
# V - the matrix containing the drift vector evaluated at the internal stages v(Q_i)
# B - the array containing the diffusion matrix evaluated at the internal stages B(Q_i)
# bdrift - the Runge-Kutta coefficients for the drift part
# bdiff - the Runge-Kutta coefficients for the ΔW terms of the diffusion part
# bdiff2- the Runge-Kutta coefficients for the ΔZ terms of the diffusion part
# Δt - the time step
# ΔW - the increments of the Brownian motion
# ΔZ - the integrals of the increments of the Brownian motion
function update_solution!(x::SolutionVector{T}, V::Vector{Vector{T}}, B::Vector{Matrix{T}},
                            bdrift::Vector{T}, bdiff::Vector{T}, bdiff2::Vector{T},
                            Δt::T, ΔW::Vector{T}, ΔZ::Vector{T}, Δy::Vector{T}=zero(ΔW)) where {T}

    @assert length(bdrift) == length(bdiff) == length(bdiff2) == length(V) == length(B)

    for i in eachindex(V, B)
        @assert length(x) == length(V[i]) == size(B[i], 1)
        @assert length(ΔW) == length(ΔZ) == size(B[i], 2)
    end

    local Δx::T

    # Contribution from the drift part
    for k in eachindex(x)
        Δx = 0
        for i in eachindex(bdrift, V)
            Δx += bdrift[i] * V[i][k]
        end
        x[k] += Δt * Δx
    end

    # Contribution from the diffusion part (ΔW terms)
    for k in eachindex(x)
        Δy .= 0
        for i in eachindex(bdiff, B)
            for l in eachindex(Δy)
                Δy[l] += bdiff[i] * B[i][k,l]
            end
        end
        x[k] += dot(Δy,ΔW)
    end

    # Contribution from the diffusion part (ΔZ terms)
    for k in eachindex(x)
        Δy .= 0
        for i in eachindex(bdiff2, B)
            for l in eachindex(Δy)
                Δy[l] += bdiff2[i] * B[i][k,l]
            end
        end
        x[k] += dot(Δy,ΔZ)/Δt
    end
end


# For stochastic partitioned Runge-Kutta methods
# q, p - the solution vector to be updated
# V, F - the matrix containing the drift vectors evaluated at the internal stages v(Q_i), f(Q_i)
# B, G - the array containing the diffusion matrices evaluated at the internal stages B(Q_i), G(Q_i)
# bqdrift, bpdrift - the Runge-Kutta coefficients for the drift parts of the q and p equations
# bqdiff, bpdiff - the Runge-Kutta coefficients for the diffusion parts of the q and p equations
# Δt - the time step
# ΔW - the increments of the Brownian motion
function update_solution!(q::SolutionVector{T}, p::SolutionVector{T},
                          V::Vector{Vector{T}}, F::Vector{Vector{T}},
                          B::Vector{Matrix{T}}, G::Vector{Matrix{T}},
                          bqdrift::Vector{T}, bqdiff::Vector{T},
                          bpdrift::Vector{T}, bpdiff::Vector{T},
                          Δt::T, ΔW::Vector{T}, Δy::Vector{T}=zero(ΔW), Δz::Vector{T}=zero(ΔW)) where {T}

    @assert length(bqdrift) == length(bqdiff) == length(bpdrift) == length(bpdiff) == length(V) == length(F) == length(B) == length(G)

    for i in eachindex(V, F, B, G)
        @assert length(q) == length(p) == length(V[i]) == length(F[i]) == size(B[i], 1) == size(G[i], 1)
        @assert length(ΔW) == size(B[i], 2) == size(G[i], 2)
    end

    local Δq::T
    local Δp::T

    # Contribution from the drift part
    for k in eachindex(q, p)
        Δq = 0
        Δp = 0
        for i in eachindex(bqdrift, bpdrift, V, F)
            Δq += bqdrift[i] * V[i][k]
            Δp += bpdrift[i] * F[i][k]
        end
        q[k] += Δt * Δq
        p[k] += Δt * Δp
    end

    # Contribution from the diffusion part
    for k in eachindex(q, p)
        Δy .= 0
        Δz .= 0

        for i in eachindex(bqdiff, bpdiff, B, G)
            for l in eachindex(Δy, Δz)
                Δy[l] += bqdiff[i] * B[i][k,l]
                Δz[l] += bpdiff[i] * G[i][k,l]
            end
        end

        q[k] += dot(Δy,ΔW)
        p[k] += dot(Δz,ΔW)
    end
end


# For stochastic split partitioned Runge-Kutta methods
# q, p - the solution vector to be updated
# V, F1, F2 - the matrix containing the drift vectors evaluated at the internal stages v(Q_i), fi(Q_i)
# B, G1, G2 - the array containing the diffusion matrices evaluated at the internal stages B(Q_i), Gi(Q_i)
# bqdrift, bpdrift1, bpdrift2 - the Runge-Kutta coefficients for the drift parts of the q and p equations
# bqdiff, bpdiff1, bpdiff2 - the Runge-Kutta coefficients for the diffusion parts of the q and p equations
# Δt - the time step
# ΔW - the increments of the Brownian motion
function update_solution!(q::SolutionVector{T}, p::SolutionVector{T},
                            V::Vector{Vector{T}}, F1::Vector{Vector{T}}, F2::Vector{Vector{T}},
                            B::Vector{Matrix{T}}, G1::Vector{Matrix{T}}, G2::Vector{Matrix{T}},
                            bqdrift::Vector{T}, bqdiff::Vector{T},
                            bpdrift1::Vector{T}, bpdrift2::Vector{T},
                            bpdiff1::Vector{T}, bpdiff2::Vector{T},
                            Δt::T, ΔW::Vector{T}, Δy::Vector{T}=zero(ΔW), Δz::Vector{T}=zero(ΔW)) where {T}

   @assert length(bqdrift) == length(bqdiff) == length(bpdrift1) == length(bpdrift2) == length(bpdiff1) == length(bpdiff2) == length(V) == length(F1) == length(F2) == length(B) == length(G1) == length(G2)

   for i in eachindex(V, F1, F2, B, G1, G2)
       @assert length(q) == length(p) == length(V[i]) == length(F1[i]) == length(F2[i]) == size(B[i], 1) == size(G1[i], 1) == size(G2[i], 1)
       @assert length(ΔW) == size(B[i], 2) == size(G1[i], 2) == size(G2[i], 2)
   end

   local Δq::T
   local Δp::T

   # Contribution from the drift part
   for k in eachindex(q, p)
       Δq = 0
       Δp = 0
       for i in eachindex(bqdrift, bpdrift1, bpdrift2, V, F1, F2)
           Δq += bqdrift[i] * V[i][k]
           Δp += bpdrift1[i] * F1[i][k] + bpdrift2[i] * F2[i][k]
       end
       q[k] += Δt * Δq
       p[k] += Δt * Δp
   end

   # Contribution from the diffusion part
   for k in eachindex(q, p)
       Δy .= 0
       Δz .= 0

       for i in eachindex(bqdiff, bpdiff1, bpdiff2, B, G1, G2)
           for l in eachindex(Δy,Δz)
               Δy[l] += bqdiff[i]  * B[i][k,l]
               Δz[l] += bpdiff1[i] * G1[i][k,l] + bpdiff2[i] * G2[i][k,l]
           end
       end

       q[k] += dot(Δy,ΔW)
       p[k] += dot(Δz,ΔW)
   end
end


# For weak Runge-Kutta methods WERK
# x - the solution vector to be updated
# V - the matrix containing the drift vector evaluated at the internal stages v(Q_i)
# B1 - the array containing the diffusion matrix evaluated at the internal stages H^(l)_i, such that B1[:,l,i] is evaluated at H^(l)_i
# B2 - the array containing the diffusion matrix evaluated at the internal stages \hat H^(l)_i, such that B2[:,l,i] is evaluated at \hat H^(l)_i
# α - the Runge-Kutta coefficients for the drift part
# β1 - the Runge-Kutta coefficients for the diffusion term with the random increments
# β2 - the Runge-Kutta coefficients for the second diffusion term
# Δt - the time step
# ΔW - the increments of the Brownian motion represented by the random variables \hat I^(k)
function update_solution!(x::SolutionVector{T}, V::Vector{Vector{T}},
                            B1::Vector{Matrix{T}}, B2::Vector{Matrix{T}}, α::Vector{T},
                            β1::Vector{T}, β2::Vector{T}, Δt::T, ΔW::Vector{T}, Δy::Vector{T}=zero(ΔW)) where {T}
    @assert length(α) == length(β1) == length(β2) == length(V) == length(B1) == length(B2)

    for i in eachindex(V, B1, B2)
        @assert length(x) == length(V[i]) == size(B1[i], 1) == size(B2[i], 1)
        @assert length(ΔW)== size(B1[i], 2) == size(B2[i], 2)
    end

    local Δx::T

    # Contribution from the drift part
    for k in eachindex(x)
        Δx = 0
        for i in eachindex(α, V)
            Δx += α[i] * V[i][k]
        end
        x[k] += Δt * Δx
    end

    # Contribution from the diffusion term with the random variables I^(k)_i
    for k in eachindex(x)
        Δy .= 0
        for i in eachindex(β1, B1)
            for l in eachindex(Δy)
                Δy[l] += β1[i] * B1[i][k,l]
            end
        end
        x[k] += dot(Δy,ΔW)
    end

    # Contribution from the second diffusion term
    for k in eachindex(x)
        Δx = 0
        for i in eachindex(β2, B2)
            for l in axes(B2[i], 2)
                Δx += β2[i] * B2[i][k,l]
            end
        end
        x[k] += sqrt(Δt)*Δx
    end
end


# For stochastic Runge-Kutta methods
function update_solution!(x::SolutionVector{T}, V::Vector{Vector{T}}, B::Vector{Matrix{T}},
                            bdrift::Vector{T}, b̂drift::Vector, bdiff::Vector{T}, b̂diff::Vector,
                            Δt::T, ΔW::Vector{T}, Δy::Vector{T}=zero(ΔW)) where {T}
    update_solution!(x, V, B, bdrift, bdiff, Δt, ΔW, Δy)
    update_solution!(x, V, B, b̂drift, b̂diff, Δt, ΔW, Δy)
end

# For stochastic partitioned Runge-Kutta methods
function update_solution!(q::SolutionVector{T}, p::SolutionVector{T},
                            V::Vector{Vector{T}}, F::Vector{Vector{T}}, B::Vector{Matrix{T}}, G::Vector{Matrix{T}},
                            bqdrift::Vector{T}, b̂qdrift::Vector, bqdiff::Vector{T}, b̂qdiff::Vector,
                            bpdrift::Vector{T}, b̂pdrift::Vector, bpdiff::Vector{T}, b̂pdiff::Vector,
                            Δt::T, ΔW::Vector{T}, Δy::Vector{T}=zero(ΔW), Δz::Vector{T}=zero(ΔW)) where {T}
    update_solution!(q, p, V, F, B, G, bqdrift, bqdiff, bpdrift, bpdiff, Δt, ΔW, Δy, Δz)
    update_solution!(q, p, V, F, B, G, b̂qdrift, b̂qdiff, b̂pdrift, b̂pdiff, Δt, ΔW, Δy, Δz)
end

 # For weak Runge-Kutta methods WERK
 function update_solution!(x::SolutionVector{T}, V::Vector{Vector{T}}, B1::Vector{Matrix{T}}, B2::Vector{Matrix{T}},
                             α::Vector{T}, α̂::Vector, β1::Vector{T}, β̂1::Vector, β2::Vector{T}, β̂2::Vector,
                             Δt::T, ΔW::Vector{T}, Δy::Vector{T}=zero(ΔW)) where {T}
     update_solution!(x, V, B1, B2, α, β1, β2, Δt, ΔW, Δy)
     update_solution!(x, V, B1, B2, α̂, β̂1, β̂2, Δt, ΔW, Δy)
 end
