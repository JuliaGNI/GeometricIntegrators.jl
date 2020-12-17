
function rel_err(sol, ref)
    maximum(abs.((sol.d[end,begin] .- ref) ./ ref))
end

function l2_err(sol, ref)
    sqrt(sum((sol.d[end,begin] .- ref).^2) / sum(ref.^2))
end
