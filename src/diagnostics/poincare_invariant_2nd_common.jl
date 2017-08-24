
using ApproxFun


@generated function compute_integral(Is)
    SU  = Ultraspherical(1, 0..1)^2
    Q   = DefiniteIntegral(SU[1]) ⊗ DefiniteIntegral(SU[2])
    # Q   = DefiniteIntegral(SU)

    quote
        return Number( $Q * Fun($SU, ApproxFun.transform($SU, Is)) )
    end
end

@generated function compute_derivative1(q)
    Dx = Derivative(Chebyshev(0..1)^2, [1,0])
    CU = Conversion(Ultraspherical(1, 0..1) ⊗ Chebyshev(0..1), Ultraspherical(1, 0..1)^2)

    quote
        return values($CU * ($Dx * q))
    end
end

@generated function compute_derivative2(q)
    Dy = Derivative(Chebyshev(0..1)^2, [0,1])
    CU = Conversion(Chebyshev(0..1) ⊗ Ultraspherical(1, 0..1), Ultraspherical(1, 0..1)^2)

    quote
        return values($CU * ($Dy * q))
    end
end

@generated function compute_canonical_invariant(pinv::Union{PoincareInvariant2nd{DT,ND,NC,NV},PoincareInvariant2ndCanonical{DT,ND,NC,NV}}, q::AbstractArray{DT,2}, p::AbstractArray{DT,2}) where {DT,ND,NC,NV}
    SC = Chebyshev(0..1)^2
    SU = Ultraspherical(1, 0..1)^2

    Js = zeros(DT, NV)
    qσ = zeros(DT, ND, NV)
    qτ = zeros(DT, ND, NV)
    pσ = zeros(DT, ND, NV)
    pτ = zeros(DT, ND, NV)

    quote
        for j in 1:size(p,1)
            fq = Fun($SC, ApproxFun.transform($SC, q[j,:]))
            fp = Fun($SC, ApproxFun.transform($SC, p[j,:]))

            # compute derivatives of q and store values
            $qσ[j,:] .= compute_derivative1(fq)
            $qτ[j,:] .= compute_derivative2(fq)

            # compute derivatives of p and store values
            $pσ[j,:] .= compute_derivative1(fp)
            $pτ[j,:] .= compute_derivative2(fp)
        end

        # compute integrands of integral invariants at all points
        for k in 1:NV
            $Js[k] = dot($pσ[:,k], $qτ[:,k]) - dot($qσ[:,k], $pτ[:,k])
        end

        # compute canonical integral invariant
        return compute_integral($Js)
    end
end

@generated function compute_noncanonical_invariant(pinv::PoincareInvariant2nd{DT,ND,NC,NV}, t::DT, q::AbstractArray{DT,2}) where {DT,ND,NC,NV}
    SC = Chebyshev(0..1)^2
    SU = Ultraspherical(1, 0..1)^2
    CU = Conversion(SC, SU)

    Is = zeros(DT, NV)
    qs = zeros(DT, ND, NV)
    qσ = zeros(DT, ND, NV)
    qτ = zeros(DT, ND, NV)
    Ω  = zeros(DT, ND, ND)

    quote
        for j in 1:size(q,1)
            fq = Fun($SC, ApproxFun.transform($SC, q[j,:]))

            # obtain values of q at the same points as the derivatives
            $qs[j,:] .= itransform($SU, resize!(($CU * fq).coefficients, NC))

            # compute derivatives of q and store values
            $qσ[j,:] .= compute_derivative1(fq)
            $qτ[j,:] .= compute_derivative2(fq)
        end

        # compute integrands of integral invariants at all points
        for k in 1:NV
            pinv.ω(t, $qs[:,k], $Ω)
            $Is[k] = vector_matrix_vector_product($qτ[:,k], $Ω, $qσ[:,k])
        end

        # compute noncanonical integral invariant
        return compute_integral($Is)
    end
end


@generated function compute_noncanonical_correction(pinv::PoincareInvariant2nd{DT,ND,NC,NV}, t::DT, q::AbstractArray{DT,2}, λ::AbstractArray{DT,2}) where {DT,ND,NC,NV}
    SC = Chebyshev(0..1)^2
    SU = Ultraspherical(1, 0..1)^2
    CU = Conversion(SC, SU)

    Ks = zeros(DT, NV)
    qs = zeros(DT, ND, NV)
    λs = zeros(DT, ND, NV)
    qσ = zeros(DT, ND, NV)
    qτ = zeros(DT, ND, NV)
    λσ = zeros(DT, ND, NV)
    λτ = zeros(DT, ND, NV)
    Ω  = zeros(DT, ND, ND)
    D²ϑ= zeros(DT, ND, ND)

    quote
        for j in 1:size(q,1)
            fq = Fun($SC, ApproxFun.transform($SC, q[j,:]))
            fλ = Fun($SC, ApproxFun.transform($SC, λ[j,:]))

            # obtain values of q at the same points as the derivatives
            $qs[j,:] .= itransform($SU, resize!(($CU * fq).coefficients, NC))
            $λs[j,:] .= itransform($SU, resize!(($CU * fλ).coefficients, NC))

            # compute derivatives of q and store values
            $qσ[j,:] .= compute_derivative1(fq)
            $qτ[j,:] .= compute_derivative2(fq)

            # compute derivatives of λ and store values
            $λσ[j,:] .= compute_derivative1(fλ)
            $λτ[j,:] .= compute_derivative2(fλ)
        end

        # compute integrands of integral invariants at all points
        for k in 1:NV
            pinv.ω(t, $qs[:,k], $Ω)
            $Ks[k] = vector_matrix_vector_product($λτ[:,k], $Ω, $λσ[:,k])

            for j in 1:length(pinv.D²ϑ)
                pinv.D²ϑ[j](t, $qs[:,k], $D²ϑ)
                $Ks[k] += 2 * $λs[j,k] * vector_matrix_vector_product($qτ[:,k], $D²ϑ, $λσ[:,k])
                $Ks[k] -= 2 * $λs[j,k] * vector_matrix_vector_product($qσ[:,k], $D²ϑ, $λτ[:,k])
            end
        end

        # compute correction to noncanonical integral invariant
        return compute_integral($Ks)
    end
end
