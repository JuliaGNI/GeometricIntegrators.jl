# Starts that differ from `q₀` by k = 0, …, K ULPs in one component, each component in turn, so
# `length(q₀) * (K + 1)` starts. A result that depends on rounding is asserted as a rate over these
# starts, not at `q₀` alone: a CPU with another SIMD width rounds like a start a few ULPs away.
function perturbed_starts(q₀, K)
    return [[i == c ? nextfloat(q₀[i], k) : q₀[i] for i in eachindex(q₀)]
            for c in eachindex(q₀) for k in 0:K]
end
