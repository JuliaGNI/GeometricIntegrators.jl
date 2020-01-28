
function rel_err(sol, ref)
    maximum(abs.((sol.d[:,end] .- ref) ./ ref))
end

function l2_err(sol, ref)
    sqrt(sum((sol.d[:,end] .- ref).^2) / sum(ref.^2))
end
