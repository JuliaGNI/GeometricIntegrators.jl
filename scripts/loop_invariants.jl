# Shared loop machinery for the Poincaré-invariant verification scripts.
#
# `include`d by `slrk_verification.jl` and `vspark_projection_symplecticity.jl`.
#
# The first Poincaré invariant of a map on (q,p) is the loop integral ∮p·dq. For any
# symplectic map it is exactly conserved, whatever the constraint structure — no
# constraint manifold and no finite-difference Jacobian are needed, which makes it the
# cleanest available symplecticity test.
#
# dq/ds along the loop is differentiated *spectrally*. A central-difference stencil
# leaves an O(N⁻²) quadrature error of order 1e-7 that masquerades as a symplecticity
# defect; the spectral derivative of an analytic closed loop converges geometrically in
# N, so the quadrature sits at round-off and the measured drift is the method's.

"""
    fourier_diffmatrix(N)

Fourier spectral differentiation matrix for `N` uniform points on `[0,2π)`.

Only the even-`N` (`cot`) form is implemented — for odd `N` the entries are built from
`csc` instead, and using this matrix would give a silently wrong derivative.
"""
function fourier_diffmatrix(N)
    iseven(N) || throw(ArgumentError("fourier_diffmatrix is the even-N (cot) form, got N = $N"))
    h = 2π / N
    D = zeros(N, N)
    for j in 1:N, k in 1:N
        j == k && continue
        D[j, k] = 0.5 * (-1)^(j + k) / tan((j - k) * h / 2)
    end
    D
end

"""
    loop_integral(qs, ps, D = fourier_diffmatrix(length(qs)))

∮ p·dq on a closed loop sampled at `N` uniform values of the parameter `s ∈ [0,1)`,
with `dq/ds` taken spectrally. For a smooth loop this is accurate to round-off.
"""
function loop_integral(qs, ps, D = fourier_diffmatrix(length(qs)))
    N = length(qs)
    acc = 0.0
    for c in eachindex(qs[1])
        dqc = 2π .* (D * [q[c] for q in qs])
        acc += sum(ps[i][c] * dqc[i] for i in 1:N) / N
    end
    acc
end

"loop of `n` points, a circle of radius `r` in the (q₁,q₂) plane through `c`"
function circle(c, r, n)
    map(1:n) do i
        q = collect(float.(c))
        q[1] += r * cos(2π * (i - 1) / n)
        q[2] += r * sin(2π * (i - 1) / n)
        q
    end
end
