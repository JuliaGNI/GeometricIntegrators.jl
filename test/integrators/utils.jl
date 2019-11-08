

function rel_energy_err_sde(sol)

    NQ = ndims(sol.q)

    if NQ==2

        en_ref  = 0.5 * ( sol.q.d[1,1]^2 + sol.q.d[2,1]^2 )
        en_last = 0.5 * ( sol.q.d[1,end]^2 + sol.q.d[2,end]^2 )

        return abs((en_last-en_ref)/en_ref)

    elseif NQ==3

        en_ref  = 0.5 * ( sol.q.d[1,1,:].^2 + sol.q.d[2,1,:].^2 )
        en_last = 0.5 * ( sol.q.d[1,end,:].^2 + sol.q.d[2,end,:].^2 )

        return maximum(abs.( (en_last .- en_ref) ./ en_ref ))

    elseif NQ==4

        en_ref  = 0.5 * ( sol.q.d[1,1,:,:].^2 + sol.q.d[2,1,:,:].^2 )
        en_last = 0.5 * ( sol.q.d[1,end,:,:].^2 + sol.q.d[2,end,:,:].^2 )

        return maximum(abs.( (en_last .- en_ref) ./ en_ref ))

    end
end


function rel_energy_err_psde(sol)

    NQ = ndims(sol.q)

    if NQ==2

        en_ref  = 0.5 * ( sol.q.d[1,1]^2 + sol.p.d[1,1]^2 )
        en_last = 0.5 * ( sol.q.d[1,end]^2 + sol.p.d[1,end]^2 )

        return abs((en_last-en_ref)/en_ref)

    elseif NQ==3

        en_ref  = 0.5 * ( sol.q.d[1,1,:].^2 + sol.p.d[1,1,:].^2 )
        en_last = 0.5 * ( sol.q.d[1,end,:].^2 + sol.p.d[1,end,:].^2 )

        return maximum(abs.( (en_last .- en_ref) ./ en_ref ))

    elseif NQ==4

        en_ref  = 0.5 * ( sol.q.d[1,1,:,:].^2 + sol.p.d[1,1,:,:].^2 )
        en_last = 0.5 * ( sol.q.d[1,end,:,:].^2 + sol.p.d[1,end,:,:].^2 )

        return maximum(abs.( (en_last .- en_ref) ./ en_ref ))

    end
end
