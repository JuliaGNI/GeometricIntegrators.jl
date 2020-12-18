using GeometricIntegrators.Utils

function rel_energy_err(sol::SolutionSDE{AT,TT,WT,1}) where {AT,TT,WT}
    en_ref  = 0.5 * L2norm(sol.q[begin])
    en_last = 0.5 * L2norm(sol.q[end])

    return abs((en_last-en_ref)/en_ref)
end

function rel_energy_err(sol::SolutionSDE{AT,TT,WT,2}) where {AT,TT,WT}
    en_ref  = 0.5 .* map(q -> L2norm(q), sol.q[begin,:])
    en_last = 0.5 .* map(q -> L2norm(q), sol.q[end,:])

    return maximum(abs.( (en_last .- en_ref) ./ en_ref ))
end

function rel_energy_err(sol::SolutionPSDE{AT,TT,WT,1}) where {AT,TT,WT}
    en_ref  = 0.5 * ( L2norm(sol.q[begin]) + L2norm(sol.p[begin]) )
    en_last = 0.5 * ( L2norm(sol.q[end])   + L2norm(sol.p[end])   )

    return abs((en_last-en_ref)/en_ref)
end

function rel_energy_err(sol::SolutionPSDE{AT,TT,WT,2}) where {AT,TT,WT}
    en_ref  = 0.5 .* map(q -> L2norm(q), sol.q[begin,:]) .+
              0.5 .* map(p -> L2norm(p), sol.p[begin,:])
    en_last = 0.5 .* map(q -> L2norm(q), sol.q[end,:]) .+
              0.5 .* map(p -> L2norm(p), sol.p[end,:])

    return maximum(abs.( (en_last .- en_ref) ./ en_ref ))
end
