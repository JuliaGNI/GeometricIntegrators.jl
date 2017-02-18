
const R₀ = 1
const B₀ = 1
const q  = 2


function dA3d3(t,x)
    zero(eltype(x))
end

function dA3d1(t,x)
    0.5*B₀*(R₀^2 - x[1]^2 + x[2]^2)/(q*x[1]^2)
end

function b3(t,x)
    -B₀*R₀*q/(sqrt(R₀^2*q^2 + x[2]^2 + (R₀ - x[1])^2)*B₀)
end

function B3(t,x)
    -B₀*R₀/x[1]
end

function dB2d3(t,x)
    zero(eltype(x))
end

function dB2d1(t,x)
    B₀/(q*x[1]) - B₀*(-R₀ + x[1])/(q*x[1]^2)
end

function db3d3(t,x)
    zero(eltype(x))
end

function db1d1(t,x)
    -B₀*x[2]*(R₀ - x[1])*q/(q*(R₀^2*q^2 + x[2]^2 + (R₀ - x[1])^2)^(3/2)*B₀)
end

function dA1d2(t,x)
    0.5*B₀*R₀/x[1]
end

function dA2d3(t,x)
    zero(eltype(x))
end

function dA1d3(t,x)
    zero(eltype(x))
end

function B1(t,x)
    -B₀*x[2]/(q*x[1])
end

function A3(t,x)
    -0.5*B₀*(x[2]^2 + (-R₀ + x[1])^2)/(q*x[1])
end

function dB3d3(t,x)
    zero(eltype(x))
end

function dB1d2(t,x)
    -B₀/(q*x[1])
end

function db1d2(t,x)
    B₀*x[2]^2*q/(q*(R₀^2*q^2 + x[2]^2 + (R₀ - x[1])^2)^(3/2)*B₀) - B₀*q/(q*sqrt(R₀^2*q^2 + x[2]^2 + (R₀ - x[1])^2)*B₀)
end

function A2(t,x)
    -0.5*B₀*R₀*log(x[1]/R₀)
end

function dBd1(t,x)
    (-R₀ + x[1])*B₀/(x[1]*sqrt(R₀^2*q^2 + x[2]^2 + (R₀ - x[1])^2)*q) - sqrt(R₀^2*q^2 + x[2]^2 + (R₀ - x[1])^2)*B₀/(x[1]^2*q)
end

function B2(t,x)
    B₀*(-R₀ + x[1])/(q*x[1])
end

function phi(t,x)
    x[3]
end

function dA2d1(t,x)
    -0.5*B₀*R₀/x[1]
end

function dB1d3(t,x)
    zero(eltype(x))
end

function db2d1(t,x)
    -B₀*(R₀ - x[1])^2*q/(q*(R₀^2*q^2 + x[2]^2 + (R₀ - x[1])^2)^(3/2)*B₀) + B₀*q/(q*sqrt(R₀^2*q^2 + x[2]^2 + (R₀ - x[1])^2)*B₀)
end

function dA3d2(t,x)
    -B₀*x[2]/(q*x[1])
end

function b1(t,x)
    -B₀*x[2]*q/(q*sqrt(R₀^2*q^2 + x[2]^2 + (R₀ - x[1])^2)*B₀)
end

function B(t,x)
    sqrt(R₀^2*q^2 + x[2]^2 + (R₀ - x[1])^2)*B₀/(x[1]*q)
end

function db3d1(t,x)
    -B₀*R₀*(R₀ - x[1])*q/((R₀^2*q^2 + x[2]^2 + (R₀ - x[1])^2)^(3/2)*B₀)
end

function dB3d2(t,x)
    zero(eltype(x))
end

function dBd3(t,x)
    zero(eltype(x))
end

function dB2d2(t,x)
    zero(eltype(x))
end

function R(t,x)
    x[1]
end

function db3d2(t,x)
    B₀*R₀*x[2]*q/((R₀^2*q^2 + x[2]^2 + (R₀ - x[1])^2)^(3/2)*B₀)
end

function Z(t,x)
    x[2]
end

function dB1d1(t,x)
    B₀*x[2]/(q*x[1]^2)
end

function b2(t,x)
    -B₀*(R₀ - x[1])*q/(q*sqrt(R₀^2*q^2 + x[2]^2 + (R₀ - x[1])^2)*B₀)
end

function dBd2(t,x)
    x[2]*B₀/(x[1]*sqrt(R₀^2*q^2 + x[2]^2 + (R₀ - x[1])^2)*q)
end

function dB3d1(t,x)
    B₀*R₀/x[1]^2
end

function db2d2(t,x)
    B₀*x[2]*(R₀ - x[1])*q/(q*(R₀^2*q^2 + x[2]^2 + (R₀ - x[1])^2)^(3/2)*B₀)
end

function dA1d1(t,x)
    -0.5*B₀*R₀*x[2]/x[1]^2
end

function db1d3(t,x)
    zero(eltype(x))
end

function A1(t,x)
    0.5*B₀*R₀*x[2]/x[1]
end

function dA2d2(t,x)
    zero(eltype(x))
end

function db2d3(t,x)
    zero(eltype(x))
end
