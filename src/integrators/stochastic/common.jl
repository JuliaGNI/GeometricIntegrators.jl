
function truncate_increments!(ΔW, A)
    if A > 0
        for i in eachindex(ΔW)
            if ΔW[i] < -A
                ΔW[i] = -A
            elseif ΔW[i] > A
                ΔW[i] = A
            end
        end
    end
end

function update_params!(int::IntegratorSIRK, sol::AtomicSolutionSDE)
    # set time for nonlinear solver and copy previous solution
    int.params.t  = sol.t
    int.params.q .= sol.q
    int.params.ΔW .= sol.ΔW
    int.params.ΔZ .= sol.ΔZ

    # truncate the increments ΔW with A
    truncate_increments!(int.params.ΔW, int.params.A)
end

function update_params!(int::Union{IntegratorSIPRK,IntegratorSISPRK}, sol::AtomicSolutionPSDE)
    # set time for nonlinear solver and copy previous solution
    int.params.t  = sol.t
    int.params.q .= sol.q
    int.params.p .= sol.p
    int.params.ΔW .= sol.ΔW
    int.params.ΔZ .= sol.ΔZ

    # truncate the increments ΔW with A
    truncate_increments!(int.params.ΔW, int.params.A)
end

function update_params!(int::IntegratorWIRK, sol::AtomicSolutionSDE)
    # set time for nonlinear solver and copy previous solution
    int.params.t  = sol.t
    int.params.q .= sol.q
    int.params.ΔW .= sol.ΔW
    int.params.ΔZ .= sol.ΔZ
end


@doc raw"""
Update solution for stochastic Runge-Kutta methods (SIRK and WIRK)

- `x`: the solution vector to be updated
- `V`: the matrix containing the drift vector evaluated at the internal stages v(Q_i) (SIRK) or v(Q0_i) (WIRK)
- `B`: the array containing the diffusion matrix evaluated at the internal stages B(Q_i) (SIRK) or B(Q1^(l)_i) (WIRK)
- `bdrift`: the Runge-Kutta coefficients for the drift part
- `bdiff`:  the Runge-Kutta coefficients for the diffusion part
- `Δt`: the time step
- `ΔW`: the increments of the Brownian motion (SFIRK) or the increments represented by the random variables Î^(k) (WFIRK)
"""
function update_solution!(sol::AtomicSolutionSDE{T}, V::Vector{Vector{T}}, B::Vector{Matrix{T}},
                            bdrift::Vector{T}, bdiff::Vector{T}, Δt::T, ΔW::Vector{T}, Δy::Vector{T}=zero(ΔW)) where {T}

    @assert length(bdrift) == length(bdiff) == length(V) == length(B)

    for i in eachindex(V, B)
        @assert length(sol.q) == length(V[i]) == size(B[i], 1)
        @assert length(ΔW)== size(B[i], 2)
    end

    local Δx::T

    # Contribution from the drift part
    for k in eachindex(sol.q)
        Δx = 0
        for i in eachindex(bdrift, V)
            Δx += bdrift[i] * V[i][k]
        end
        update!(sol, Δt*Δx, k)
    end

    # Contribution from the diffusion part
    for k in eachindex(sol.q)
        Δy .= 0
        for i in eachindex(bdiff, B)
            for l in eachindex(Δy)
                Δy[l] += bdiff[i] * B[i][k,l]
            end
        end
        update!(sol, dot(Δy,ΔW), k)
    end
end


@doc raw"""
Update solution for stochastic Runge-Kutta methods (SERK)

- `x`: the solution vector to be updated
- `V`: the matrix containing the drift vector evaluated at the internal stages v(Q_i)
- `B`: the array containing the diffusion matrix evaluated at the internal stages B(Q_i)
- `bdrift`: the Runge-Kutta coefficients for the drift part
- `bdiff`:  the Runge-Kutta coefficients for the ΔW terms of the diffusion part
- `bdiff2`: the Runge-Kutta coefficients for the ΔZ terms of the diffusion part
- `Δt`: the time step
- `ΔW`: the increments of the Brownian motion
- `ΔZ`: the integrals of the increments of the Brownian motion
"""
function update_solution!(sol::AtomicSolutionSDE{T}, V::Vector{Vector{T}}, B::Vector{Matrix{T}},
                            bdrift::Vector{T}, bdiff::Vector{T}, bdiff2::Vector{T},
                            Δt::T, ΔW::Vector{T}, ΔZ::Vector{T}, Δy::Vector{T}=zero(ΔW)) where {T}

    @assert length(bdiff2) == length(B)

    for i in eachindex(B)
        @assert length(ΔZ)== size(B[i], 2)
    end

    # Contributions from the drift and diffusion parts
    update_solution!(sol, V, B, bdrift, bdiff, Δt, ΔW, Δy)

    # Contribution from the diffusion part (ΔZ terms)
    for k in eachindex(sol.q)
        Δy .= 0
        for i in eachindex(bdiff2, B)
            for l in eachindex(Δy)
                Δy[l] += bdiff2[i] * B[i][k,l]
            end
        end
        update!(sol, dot(Δy,ΔZ)/Δt, k)
    end
end


@doc raw"""
Update solution for stochastic partitioned Runge-Kutta methods

- `q`, `p`: the solution vector to be updated
- `V`, `F`: the matrix containing the drift vectors evaluated at the internal stages v(Q_i), f(Q_i)
- `B`, `G`: the array containing the diffusion matrices evaluated at the internal stages B(Q_i), G(Q_i)
- `bqdrift`, `bpdrift`: the Runge-Kutta coefficients for the drift parts of the q and p equations
- `bqdiff`, `bpdiff`:   the Runge-Kutta coefficients for the diffusion parts of the q and p equations
- `Δt`: the time step
- `ΔW`: the increments of the Brownian motion
"""
function update_solution!(sol::AtomicSolutionPSDE{T},
                          V::Vector{Vector{T}}, F::Vector{Vector{T}},
                          B::Vector{Matrix{T}}, G::Vector{Matrix{T}},
                          bqdrift::Vector{T}, bqdiff::Vector{T},
                          bpdrift::Vector{T}, bpdiff::Vector{T},
                          Δt::T, ΔW::Vector{T}, Δy::Vector{T}=zero(ΔW), Δz::Vector{T}=zero(ΔW)) where {T}

    @assert length(bqdrift) == length(bqdiff) == length(bpdrift) == length(bpdiff) == length(V) == length(F) == length(B) == length(G)

    for i in eachindex(V, F, B, G)
        @assert length(sol.q) == length(sol.p) == length(V[i]) == length(F[i]) == size(B[i], 1) == size(G[i], 1)
        @assert length(ΔW) == size(B[i], 2) == size(G[i], 2)
    end

    local Δq::T
    local Δp::T

    # Contribution from the drift part
    for k in eachindex(sol.q, sol.p)
        Δq = 0
        Δp = 0
        for i in eachindex(bqdrift, bpdrift, V, F)
            Δq += bqdrift[i] * V[i][k]
            Δp += bpdrift[i] * F[i][k]
        end
        update!(sol, Δt*Δq, Δt*Δp, k)
    end

    # Contribution from the diffusion part
    for k in eachindex(sol.q, sol.p)
        Δy .= 0
        Δz .= 0
        for i in eachindex(bqdiff, bpdiff, B, G)
            for l in eachindex(Δy, Δz)
                Δy[l] += bqdiff[i] * B[i][k,l]
                Δz[l] += bpdiff[i] * G[i][k,l]
            end
        end
        update!(sol, dot(Δy,ΔW), dot(Δz,ΔW), k)
    end
end


@doc raw"""
Update solution for stochastic split partitioned Runge-Kutta methods

- `q`, `p`: the solution vector to be updated
- `V`, `F1`, `F2`: the matrix containing the drift vectors evaluated at the internal stages v(Q_i), fi(Q_i)
- `B`, `G1`, `G2`: the array containing the diffusion matrices evaluated at the internal stages B(Q_i), Gi(Q_i)
- `bqdrift`, `bpdrift1`, `bpdrift2`: the Runge-Kutta coefficients for the drift parts of the q and p equations
- `bqdiff`, `bpdiff1`, `bpdiff2`:    the Runge-Kutta coefficients for the diffusion parts of the q and p equations
- `Δt`: the time step
- `ΔW`: the increments of the Brownian motion
"""
function update_solution!(sol::AtomicSolutionPSDE{T},
                            V::Vector{Vector{T}}, F1::Vector{Vector{T}}, F2::Vector{Vector{T}},
                            B::Vector{Matrix{T}}, G1::Vector{Matrix{T}}, G2::Vector{Matrix{T}},
                            bqdrift::Vector{T}, bqdiff::Vector{T},
                            bpdrift1::Vector{T}, bpdrift2::Vector{T},
                            bpdiff1::Vector{T}, bpdiff2::Vector{T},
                            Δt::T, ΔW::Vector{T}, Δy::Vector{T}=zero(ΔW), Δz::Vector{T}=zero(ΔW)) where {T}

   @assert length(bqdrift) == length(bqdiff) == length(bpdrift1) == length(bpdrift2) == length(bpdiff1) == length(bpdiff2) == length(V) == length(F1) == length(F2) == length(B) == length(G1) == length(G2)

   for i in eachindex(V, F1, F2, B, G1, G2)
       @assert length(sol.q) == length(sol.p) == length(V[i]) == length(F1[i]) == length(F2[i]) == size(B[i], 1) == size(G1[i], 1) == size(G2[i], 1)
       @assert length(ΔW) == size(B[i], 2) == size(G1[i], 2) == size(G2[i], 2)
   end

   local Δq::T
   local Δp::T

   # Contribution from the drift part
   for k in eachindex(sol.q, sol.p)
       Δq = 0
       Δp = 0
       for i in eachindex(bqdrift, bpdrift1, bpdrift2, V, F1, F2)
           Δq += bqdrift[i] * V[i][k]
           Δp += bpdrift1[i] * F1[i][k] + bpdrift2[i] * F2[i][k]
       end
       update!(sol, Δt*Δq, Δt*Δp, k)
   end

   # Contribution from the diffusion part
   for k in eachindex(sol.q, sol.p)
       Δy .= 0
       Δz .= 0
       for i in eachindex(bqdiff, bpdiff1, bpdiff2, B, G1, G2)
           for l in eachindex(Δy,Δz)
               Δy[l] += bqdiff[i]  * B[i][k,l]
               Δz[l] += bpdiff1[i] * G1[i][k,l] + bpdiff2[i] * G2[i][k,l]
           end
       end
       update!(sol, dot(Δy,ΔW), dot(Δz,ΔW), k)
   end
end


@doc raw"""
Update solution for weak Runge-Kutta methods WERK

- `x`:  the solution vector to be updated
- `V`:  the matrix containing the drift vector evaluated at the internal stages v(Q_i)
- `B1`: the array containing the diffusion matrix evaluated at the internal stages H^(l)_i, such that B1[:,l,i] is evaluated at H^(l)_i
- `B2`: the array containing the diffusion matrix evaluated at the internal stages Ĥ^(l)_i, such that B2[:,l,i] is evaluated at Ĥ^(l)_i
- `α`:  the Runge-Kutta coefficients for the drift part
- `β1`: the Runge-Kutta coefficients for the diffusion term with the random increments
- `β2`: the Runge-Kutta coefficients for the second diffusion term
- `Δt`: the time step
- `ΔW`: the increments of the Brownian motion represented by the random variables Î^(k)
"""
function update_solution!(sol::AtomicSolutionSDE{T}, V::Vector{Vector{T}},
                            B1::Vector{Matrix{T}}, B2::Vector{Matrix{T}}, α::Vector{T},
                            β1::Vector{T}, β2::Vector{T}, Δt::T, ΔW::Vector{T}, Δy::Vector{T}=zero(ΔW)) where {T}
    @assert length(α) == length(β1) == length(β2) == length(V) == length(B1) == length(B2)

    for i in eachindex(V, B1, B2)
        @assert length(sol.q) == length(V[i]) == size(B1[i], 1) == size(B2[i], 1)
        @assert length(ΔW)== size(B1[i], 2) == size(B2[i], 2)
    end

    local Δx::T

    # Contribution from the drift part
    for k in eachindex(sol.q)
        Δx = 0
        for i in eachindex(α, V)
            Δx += α[i] * V[i][k]
        end
        update!(sol, Δt*Δx, k)
    end

    # Contribution from the diffusion term with the random variables I^(k)_i
    for k in eachindex(sol.q)
        Δy .= 0
        for i in eachindex(β1, B1)
            for l in eachindex(Δy)
                Δy[l] += β1[i] * B1[i][k,l]
            end
        end
        update!(sol, dot(Δy,ΔW), k)
    end

    # Contribution from the second diffusion term
    for k in eachindex(sol.q)
        Δx = 0
        for i in eachindex(β2, B2)
            for l in axes(B2[i], 2)
                Δx += β2[i] * B2[i][k,l]
            end
        end
        update!(sol, sqrt(Δt)*Δx, k)
    end
end
