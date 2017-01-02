
n = 11
x = collect(linspace(0., 1., n))
y = rand(n)

V = vandermonde_matrix(x)
Vinv = vandermonde_matrix_inverse(x)

lu = LUSolver(V, y)

factorize!(lu, false)
solve!(lu)

a = *(Vinv, y)

@test maximum(abs((a .- lu.x) ./ a)) < eps(Float32)
