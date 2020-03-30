
@testset "$(rpad("Vandermonde Matrix Tests",80))" begin

    n = 11
    x = collect(range(0., stop=1., length=n))
    y = rand(n)

    V = vandermonde_matrix(x)
    Vinv = vandermonde_matrix_inverse(x)

    x = V\y
    a = *(Vinv, y)

    @test maximum(abs.((a .- x) ./ a)) < 4*eps(Float32)

end
