
function rel_err(sol, ref)
    maximum(abs.((sol.d[:,end] .- ref) ./ ref))
end
