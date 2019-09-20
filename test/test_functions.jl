
function fx(t, x, fx)
    fx .= x
end

function fq(t, q, p, fq)
    fq .= q
end

function fp(t, q, p, fp)
    fp .= q.^2
end

function fϕ(t, x, fϕ)
    fϕ .= 0
end

function gx(t, x, λ, fx)
    fx .= x
end

function gq(t, q, p, λ, fλ)
    fλ .= q
end

function gp(t, q, p, λ, fλ)
    fλ .= q.^2
end

function gϕ(t, q, p, gϕ)
    gϕ .= p - q.^2
end
