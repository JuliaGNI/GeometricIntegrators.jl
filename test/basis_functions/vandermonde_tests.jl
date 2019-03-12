
n = 11
x = collect(range(0., stop=1., length=n))
y = rand(n)

V = vandermonde_matrix(x)
Vinv = vandermonde_matrix_inverse(x)

lu = LUSolver(V, y)

factorize!(lu, false)
solve!(lu)

a = *(Vinv, y)

@test maximum(abs.((a .- lu.x) ./ a)) < eps(Float32)
