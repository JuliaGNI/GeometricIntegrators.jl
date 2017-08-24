
# using ApproxFun: Chebyshev, Ultraspherical, Conversion, points, transform
using ApproxFun
using GeometricIntegrators.Diagnostics: compute_integral, compute_derivative1, compute_derivative2


### Test parameters ###

# nx  = 20
# ny  = 20
# dtol = 1E-12
# itol = 1E-11

nx  = 100
ny  = 100
dtol = 1E-10
itol = 1E-16

scale = 1E-00
# scale = 1E-15


### Set up ApproxFun spaces ###

SC = Chebyshev(0..1)^2
SU = Ultraspherical(1, 0..1)^2
CU = Conversion(SC, SU)


### check derivatives ###

# compute x=(cos(2πτ)cos(2πσ)) on SC
pc = points(SC, nx*ny)
xc = [scale * cos(2π*pc[i][1]) * cos(2π*pc[i][2]) for i in eachindex(pc)]

# compute x=(cos(2πτ)cos(2πσ)) on SU
pu = points(SU, nx*ny)
xu = [scale * cos(2π*pu[i][1]) * cos(2π*pu[i][2]) for i in eachindex(pu)]

# create Fun and compute derivatives
fun_xc = Fun(SC, transform(SC, xc))
fun_xu = Fun(SU, transform(SU, xu))

fx1 = compute_derivative1(fun_xc)
fx2 = compute_derivative2(fun_xc)

# compute τ-derivative of x (-2πsin(2πτ)cos(2πσ))
pu = points(SU, length(fx1))
x1 = [- 2π * scale * sin(2π*pu[i][1]) * cos(2π*pu[i][2]) for i in eachindex(pu)]

# compute σ-derivative of x (-2πcos(2πτ)sin(2πσ))
pu = points(SU, length(fx2))
x2 = [- 2π * scale * cos(2π*pu[i][1]) * sin(2π*pu[i][2]) for i in eachindex(pu)]

# compare Fun derivative and analytical derivative
# println(maximum(abs.(fx1 .- x1)))
# println(maximum(abs.(fx2 .- x2)))

# the following test should fail but for some reason does not
# (with nx=ny=20 success is expected for atol=5E-13 or larger)
# @test fx1 ≈ x1 atol=1E-15
# @test fx2 ≈ x2 atol=1E-15

# this works, but is less compact
@test (fx1 .- x1) ≈ zeros(length(fx1)) atol=dtol*scale
@test (fx2 .- x2) ≈ zeros(length(fx2)) atol=dtol*scale


### check integral ###

# compute x=(cos(2πτ)cos(2πσ))
pu = points(SU, nx*ny)
x  = zeros(length(pu))

for i in eachindex(x,pu)
    x[i] .= scale * cos(2π*pu[i][1]) * cos(2π*pu[i][2])
end

@test compute_integral(x) ≈ 0.0 atol=itol*scale
